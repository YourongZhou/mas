const state = {
  taskId: null,
  tasks: [],
  snapshot: null,
  selectedTraceId: null,
  selectedFilePath: null,
  userSelectedTab: false,
  eventSource: null,
  refreshTimer: null,
  replyTaskId: null,
  isSubmitting: false,
  submitMode: null,
  submitToken: 0,
  defaultTaskText: "",
  defaultDataPath: "",
  defaultMaxTurns: "20",
};

const els = {};

document.addEventListener("DOMContentLoaded", () => {
  for (const id of [
    "taskTitle", "runtimeLine", "currentToolPill", "statusPill", "refreshButton", "newTaskButton", "taskSearch", "clearTaskSearchButton",
    "taskList", "messageStream", "scrollLatestButton", "taskForm", "taskInput", "dataPathInput", "maxTurnsInput",
    "errorLine", "detailTabs", "detailBody", "traceTab", "resultTab", "fileTab", "traceTitle", "traceMeta", "tracePlanCount", "tracePlanList", "traceInput", "traceOutput", "traceListCount", "traceList",
    "resultSummary", "resultVisuals", "resultFileList", "fileTitle", "fileMeta", "copyFilePathButton", "downloadLink", "filePreview", "planCount",
    "planList", "filesRefreshButton", "fileSearch", "clearFileSearchButton", "fileTree", "computeState", "computeGrid", "resultComputeGrid", "sendButton",
    "composerOptions"
  ]) {
    els[id] = document.getElementById(id);
  }
  state.defaultTaskText = els.taskInput.value;
  state.defaultDataPath = els.dataPathInput.value;
  state.defaultMaxTurns = els.maxTurnsInput.value;

  els.taskForm.addEventListener("submit", startTask);
  els.currentToolPill.addEventListener("click", () => {
    if (state.snapshot?.status === "needs_user_input") {
      els.taskInput.focus();
      return;
    }
    const trace = topbarTrace(state.snapshot);
    if (!trace) return;
    state.selectedTraceId = trace.id;
    markSelectedTraceRow();
    renderTraceDetail();
    activateTab("trace", { user: true });
  });
  els.refreshButton.addEventListener("click", runRefresh);
  els.filesRefreshButton.addEventListener("click", runRefresh);
  els.taskSearch.addEventListener("input", () => renderTasks(state.tasks));
  els.clearTaskSearchButton.addEventListener("click", () => {
    els.taskSearch.value = "";
    renderTasks(state.tasks);
    els.taskSearch.focus();
  });
  els.fileSearch.addEventListener("input", () => renderFiles(state.snapshot));
  els.clearFileSearchButton.addEventListener("click", () => {
    els.fileSearch.value = "";
    renderFiles(state.snapshot);
    els.fileSearch.focus();
  });
  document.querySelectorAll("[data-demo-path]").forEach((row) => {
    row.addEventListener("click", () => {
      els.dataPathInput.value = row.dataset.demoPath;
    });
    enableKeyboardActivation(row);
  });
  els.scrollLatestButton.addEventListener("click", scrollToLatest);
  els.messageStream.addEventListener("scroll", () => {
    if (shouldStickToBottom(els.messageStream)) {
      setScrollLatestVisible(false);
    }
  });
  els.newTaskButton.addEventListener("click", () => {
    history.pushState(null, "", "/");
    showNewTask();
  });
  els.detailTabs.addEventListener("click", (event) => {
    if (event.target.matches(".tab")) {
      activateTab(event.target.dataset.tab, { user: true });
    }
  });
  els.detailTabs.addEventListener("keydown", handleDetailTabKeydown);
  document.addEventListener("click", async (event) => {
    const logButton = event.target.closest("[data-open-log]");
    if (logButton) {
      await openLog();
      return;
    }
    const rerunCurrentButton = event.target.closest("[data-rerun-current]");
    if (rerunCurrentButton) {
      prepareRerunCurrent();
      return;
    }
    const failedTraceButton = event.target.closest("[data-open-failed-trace]");
    if (failedTraceButton) {
      const trace = topbarTrace(state.snapshot);
      if (!trace) return;
      state.selectedTraceId = trace.id;
      markSelectedTraceRow();
      renderTraceDetail();
      activateTab("trace", { user: true });
      return;
    }
    const button = event.target.closest("[data-copy]");
    if (button) {
      try {
        const copyValue = button.dataset.copyValue ?? document.getElementById(button.dataset.copy)?.textContent ?? "";
        await navigator.clipboard.writeText(copyValue);
        showCopyFeedback(button, "Copied");
      } catch (error) {
        showError("Unable to copy to clipboard.");
      }
    }
  });
  window.addEventListener("popstate", handleLocationChange);

  init();
});

async function init() {
  try {
    const health = await api("/api/health");
    els.runtimeLine.textContent = `${health.model || "model"} · ${health.baseUrl || "local"}`;
    await loadTasks();
    await handleLocationChange();
  } catch (error) {
    showError(error.message);
  }
}

async function handleLocationChange() {
  cancelSubmit();
  try {
    const match = location.pathname.match(/\/tasks\/([^/]+)/);
    if (match) {
      await selectTask(match[1]);
      return;
    }
    showNewTask();
  } catch (error) {
    restoreRouteAfterNavigationFailure();
    showError(error.message);
  }
}

function restoreRouteAfterNavigationFailure() {
  const visiblePath = state.taskId ? `/tasks/${state.taskId}` : "/";
  history.replaceState(null, "", visiblePath);
}

function showNewTask() {
  state.taskId = null;
  state.snapshot = null;
  state.selectedTraceId = null;
  state.selectedFilePath = null;
  state.userSelectedTab = false;
  state.replyTaskId = null;
  cancelSubmit();
  closeEvents();
  resetFilePreview();
  resetFileSearch();
  resetComposerDefaults();
  showError("");
  activateTab("trace");
  renderTasks(state.tasks);
  render();
}

async function startTask(event) {
  event.preventDefault();
  showError("");
  if (state.isSubmitting) return;
  if (state.snapshot?.status === "needs_user_input" && state.taskId) {
    await resumeTask();
    return;
  }
  if (isActiveStatus(state.snapshot?.status)) {
    showError("This task is still running. Start a new task from New if needed.");
    return;
  }
  closeEvents();
  const payload = {
    task: els.taskInput.value.trim(),
    data_path: els.dataPathInput.value.trim() || null,
    max_turns: Number(els.maxTurnsInput.value || 20),
  };
  if (!payload.task) {
    showError("Task is required.");
    return;
  }
  const submitToken = beginSubmit("create");
  try {
    const created = await api("/api/tasks", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    if (state.submitToken !== submitToken) return;
    addUserMessage(payload.task);
    history.pushState(null, "", created.redirectUrl);
    await loadTasks();
    await selectTask(created.taskId);
  } catch (error) {
    showError(error.message);
  } finally {
    finishSubmit(submitToken);
  }
}

async function resumeTask() {
  const userAnswer = els.taskInput.value.trim();
  const taskId = state.taskId;
  if (!userAnswer) {
    showError("Reply is required.");
    return;
  }
  const submitToken = beginSubmit("reply");
  try {
    await api(`/api/tasks/${taskId}/resume`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        user_answer: userAnswer,
        max_turns: Number(els.maxTurnsInput.value || 20),
      }),
    });
    if (state.submitToken !== submitToken) return;
    els.taskInput.value = "";
    await loadTasks();
    await selectTask(taskId);
  } catch (error) {
    showError(error.message);
  } finally {
    finishSubmit(submitToken);
  }
}

async function selectTask(taskId) {
  const previousTaskId = state.taskId;
  const previousSnapshot = state.snapshot;
  state.taskId = taskId;
  try {
    const snapshot = await api(`/api/tasks/${taskId}`);
    if (state.taskId !== taskId) return;
    state.snapshot = snapshot;
    state.selectedTraceId = chooseTraceId(state.snapshot, null);
    state.userSelectedTab = false;
    mergeTaskSnapshot(state.snapshot);
    syncComposerFromSnapshot(state.snapshot);
    resetFilePreview();
    resetFileSearch();
    showError("");
    activateDefaultTab(state.snapshot);
    subscribe(taskId);
    render();
    scrollToLatest();
  } catch (error) {
    state.taskId = previousTaskId;
    state.snapshot = previousSnapshot;
    render();
    throw error;
  }
}

async function refresh() {
  await loadTasks();
  const taskId = state.taskId;
  if (taskId) {
    const snapshot = await api(`/api/tasks/${taskId}`);
    if (state.taskId !== taskId) return;
    state.snapshot = snapshot;
    state.selectedTraceId = chooseTraceId(state.snapshot, state.selectedTraceId);
    mergeTaskSnapshot(state.snapshot);
  }
  render();
}

async function runRefresh() {
  showError("");
  try {
    await refresh();
  } catch (error) {
    showError(error.message);
  }
}

async function loadTasks() {
  const payload = await api("/api/tasks");
  state.tasks = payload.tasks || [];
  renderTasks(state.tasks);
}

function subscribe(taskId) {
  closeEvents();
  const source = new EventSource(`/api/tasks/${taskId}/events`);
  state.eventSource = source;
  source.addEventListener("bioagent_event", () => {
    scheduleSnapshotRefresh(taskId);
  });
  source.addEventListener("task_snapshot", (event) => {
    if (state.taskId !== taskId) return;
    state.snapshot = JSON.parse(event.data);
    state.selectedTraceId = chooseTraceId(state.snapshot, state.selectedTraceId);
    mergeTaskSnapshot(state.snapshot);
    render();
    if (isTerminalStatus(state.snapshot?.status)) {
      closeEvents();
    }
  });
  source.onerror = () => {
    if (state.taskId !== taskId) return;
    if (isTerminalStatus(state.snapshot?.status)) {
      closeEvents();
      return;
    }
    if (state.snapshot?.status === "running") {
      showError("Realtime stream disconnected. The snapshot can still be refreshed.");
    }
  };
}

function closeEvents() {
  if (state.refreshTimer) {
    clearTimeout(state.refreshTimer);
    state.refreshTimer = null;
  }
  if (state.eventSource) {
    state.eventSource.close();
    state.eventSource = null;
  }
}

function scheduleSnapshotRefresh(taskId) {
  if (state.refreshTimer) return;
  state.refreshTimer = window.setTimeout(async () => {
    state.refreshTimer = null;
    if (state.taskId !== taskId) return;
    try {
      const snapshot = await api(`/api/tasks/${taskId}`);
      if (state.taskId !== taskId) return;
      state.snapshot = snapshot;
      state.selectedTraceId = chooseTraceId(state.snapshot, state.selectedTraceId);
      mergeTaskSnapshot(state.snapshot);
      render();
      if (isTerminalStatus(state.snapshot?.status)) {
        closeEvents();
      }
    } catch (error) {
      showError(error.message);
    }
  }, 500);
}

function isTerminalStatus(status) {
  return ["completed", "failed", "needs_user_input"].includes(status);
}

function isActiveStatus(status) {
  return ["queued", "running"].includes(status);
}

function render() {
  const snapshot = state.snapshot;
  const status = snapshot?.status || "idle";
  const displayTitle = snapshot ? formatTaskDisplayTitle(snapshot) : { title: "New task" };
  activateDefaultTab(snapshot);
  els.taskTitle.textContent = displayTitle.title;
  els.taskTitle.title = snapshot?.task || snapshot?.title || "New task";
  els.statusPill.textContent = formatStatusLabel(status);
  els.statusPill.className = `pill ${statusClass(status)}`;
  els.computeState.textContent = formatStatusLabel(status);
  els.computeState.className = `pill ${statusClass(status)}`;
  renderCurrentTool(snapshot);
  renderTabCounts(snapshot);
  renderMessages(snapshot);
  renderTraces(snapshot);
  renderPlan(snapshot);
  renderFiles(snapshot);
  renderCompute(snapshot);
  renderResultSummary(snapshot);
  renderResultVisuals(snapshot);
  renderResultFileList(snapshot);
  renderComposer(snapshot);
}

function mergeTaskSnapshot(snapshot) {
  if (!snapshot?.id) return;
  const item = {
    id: snapshot.id,
    title: snapshot.title || snapshot.task || "Untitled BioAgent task",
    task: snapshot.task || "",
    status: snapshot.status || "idle",
    error: snapshot.error || "",
    dataPath: snapshot.dataPath || "",
    logPath: snapshot.compute?.logPath || "",
    createdAt: snapshot.createdAt || new Date().toISOString(),
    updatedAt: snapshot.updatedAt || new Date().toISOString(),
    messageCount: (snapshot.messages || []).length,
    traceCount: (snapshot.traces || []).length,
    resultFileCount: flattenFiles(snapshot.resultFiles || []).length,
    resultFileNames: flattenFiles(snapshot.resultFiles || []).map((file) => file.name),
  };
  const index = state.tasks.findIndex((task) => task.id === item.id);
  if (index >= 0) {
    state.tasks[index] = { ...state.tasks[index], ...item };
  } else {
    state.tasks.unshift(item);
  }
  state.tasks.sort((left, right) => new Date(right.updatedAt) - new Date(left.updatedAt));
  renderTasks(state.tasks);
}

function chooseTraceId(snapshot, currentId) {
  const traces = snapshot?.traces || [];
  if (!traces.length) return null;
  if (currentId && traces.some((trace) => trace.id === currentId)) {
    return currentId;
  }
  const running = traces.findLast((trace) => trace.status === "running");
  if (running) return running.id;
  const failed = snapshot?.status === "failed" ? traces.findLast((trace) => trace.status === "failed") : null;
  if (failed) return failed.id;
  return traces[traces.length - 1].id;
}

function renderCurrentTool(snapshot) {
  const compute = snapshot?.compute || {};
  const trace = topbarTrace(snapshot);
  const currentTool = compute.currentTool || "";
  if (snapshot?.status === "needs_user_input") {
    els.currentToolPill.textContent = "Waiting for reply";
    els.currentToolPill.className = "tool-pill needs-input";
    els.currentToolPill.disabled = false;
    els.currentToolPill.title = "Focus reply box";
    return;
  }
  const failure = formatSnapshotFailure(snapshot);
  if (failure) {
    els.currentToolPill.textContent = `Failed: ${failure}`;
    els.currentToolPill.className = "tool-pill failed";
    els.currentToolPill.disabled = !trace;
    els.currentToolPill.title = trace ? "Open failed tool trace" : "No tool trace yet";
    return;
  }
  if (!trace && !currentTool) {
    els.currentToolPill.textContent = "No tool trace";
    els.currentToolPill.className = "tool-pill idle";
    els.currentToolPill.disabled = true;
    els.currentToolPill.title = "No tool trace yet";
    return;
  }
  const status = currentTool ? "running" : trace ? formatTraceStatusClass(trace, snapshot) : "idle";
  const label = currentTool || trace?.toolName || "idle";
  const prefix = currentToolLabelPrefix(currentTool, snapshot);
  els.currentToolPill.textContent = `${prefix}: ${formatToolName(label)}`;
  els.currentToolPill.className = `tool-pill ${status}`;
  els.currentToolPill.disabled = !trace;
  els.currentToolPill.title = trace ? "Open tool trace" : "No tool trace yet";
}

function renderTabCounts(snapshot) {
  const traceCount = snapshot?.traces?.length || 0;
  const resultCount = flattenFiles(snapshot?.resultFiles || []).length;
  els.traceTab.textContent = traceCount ? `Trace ${traceCount}` : "Trace";
  els.resultTab.textContent = resultCount ? `Result ${resultCount}` : "Result";
  els.fileTab.textContent = "File";
}

function currentToolLabelPrefix(currentTool, snapshot) {
  if (currentTool) return "Running";
  if (isTerminalStatus(snapshot?.status)) return "Last tool";
  return "Tool";
}

function topbarTrace(snapshot) {
  const traces = snapshot?.traces || [];
  const running = traces.findLast((trace) => trace.status === "running");
  const failed = snapshot?.status === "failed" ? traces.findLast((trace) => trace.status === "failed") : null;
  return running || failed || traces[traces.length - 1] || null;
}

function renderTasks(tasks) {
  const query = els.taskSearch.value.trim().toLowerCase();
  syncTaskSearchClearButton(query);
  const filtered = tasks.filter((task) => taskSearchText(task).includes(query));
  els.taskList.innerHTML = filtered.length ? "" : `<div class="empty">${tasks.length ? "No matching runs." : "No runs yet."}</div>`;
  for (const task of filtered) {
    const row = document.createElement("div");
    const displayTitle = formatTaskDisplayTitle(task);
    row.className = `task-row ${task.id === state.taskId ? "active" : ""}`;
    enableKeyboardActivation(row);
    row.title = formatTaskTooltip(task);
    row.setAttribute("aria-label", formatTaskAriaLabel(task));
    row.innerHTML = `
      <div class="row-title"><strong>${escapeHtml(displayTitle.title)}</strong>${displayTitle.detail ? `<span>${escapeHtml(displayTitle.detail)}</span>` : ""}</div>
      <div class="task-actions">
        <span class="pill ${statusClass(task.status)}">${escapeHtml(formatStatusLabel(task.status))}</span>
        ${Number(task.traceCount || 0) > 0 ? `<button class="icon-button task-trace" type="button" data-open-task-trace-id="${escapeHtml(task.id)}">Trace</button>` : ""}
        ${task.logPath ? `<button class="icon-button task-log" type="button" data-open-task-log-id="${escapeHtml(task.id)}">Log</button>` : ""}
        ${isRerunnableTask(task) ? `<button class="icon-button task-rerun" type="button" data-rerun-task-id="${escapeHtml(task.id)}">Rerun</button>` : ""}
      </div>
      <div class="row-sub">${escapeHtml(formatTaskDisplayMeta(task, query))}</div>
    `;
    const traceButton = row.querySelector("[data-open-task-trace-id]");
    traceButton?.addEventListener("click", (event) => prepareOpenTraceFromList(event, task));
    const logButton = row.querySelector("[data-open-task-log-id]");
    logButton?.addEventListener("click", (event) => prepareOpenLogFromList(event, task));
    const rerunButton = row.querySelector("[data-rerun-task-id]");
    rerunButton?.addEventListener("click", (event) => prepareRerunFromList(event, task));
    row.addEventListener("click", (event) => {
      if (event.target.closest("[data-open-task-trace-id]")) return;
      if (event.target.closest("[data-open-task-log-id]")) return;
      if (event.target.closest("[data-rerun-task-id]")) return;
      if (task.id === state.taskId) return;
      cancelSubmit();
      history.pushState(null, "", `/tasks/${task.id}`);
      openTaskFromList(task.id);
    });
    els.taskList.appendChild(row);
  }
}

function syncTaskSearchClearButton(query) {
  els.clearTaskSearchButton.hidden = !query;
}

function syncFileSearchClearButton(query) {
  els.clearFileSearchButton.hidden = !query;
}

function isRerunnableTask(task) {
  return ["completed", "failed"].includes(task.status);
}

async function prepareOpenTraceFromList(event, task) {
  event.stopPropagation();
  cancelSubmit();
  if (task.id !== state.taskId) {
    history.pushState(null, "", `/tasks/${task.id}`);
    await openTaskFromList(task.id);
  }
  if (state.taskId !== task.id) return;
  activateTab("trace", { user: true });
}

async function prepareOpenLogFromList(event, task) {
  event.stopPropagation();
  cancelSubmit();
  if (task.id !== state.taskId) {
    history.pushState(null, "", `/tasks/${task.id}`);
    await openTaskFromList(task.id);
  }
  if (state.taskId !== task.id) return;
  await openLog();
}

async function prepareRerunFromList(event, task) {
  event.stopPropagation();
  cancelSubmit();
  if (task.id !== state.taskId) {
    history.pushState(null, "", `/tasks/${task.id}`);
    await openTaskFromList(task.id);
  }
  if (state.taskId !== task.id) return;
  els.taskInput.focus();
  els.taskInput.setSelectionRange?.(els.taskInput.value.length, els.taskInput.value.length);
}

function prepareRerunCurrent() {
  if (!isTerminalStatus(state.snapshot?.status)) return;
  syncComposerFromSnapshot(state.snapshot);
  renderComposer(state.snapshot);
  els.taskInput.focus();
  els.taskInput.setSelectionRange?.(els.taskInput.value.length, els.taskInput.value.length);
}

function formatTaskDisplayTitle(task) {
  const text = String(task.task || task.title || "Untitled BioAgent task").trim();
  const firstLine = text.split(/\r?\n/)[0] || text;
  const firstSentence = firstLine.split(/[。.!?]/)[0].trim() || firstLine.trim() || "Untitled BioAgent task";
  const title = firstSentence.length > 48 ? `${firstSentence.slice(0, 48)}...` : firstSentence;
  const dataName = basename(task.dataPath);
  return { title, detail: dataName };
}

function basename(path) {
  return String(path || "").split(/[\\/]/).filter(Boolean).pop() || "";
}

function formatTaskDisplayMeta(task, query = "") {
  const updated = new Date(task.updatedAt).toLocaleString();
  const duration = formatRunDuration(task);
  const activity = formatTaskActivityCounts(task);
  const outputMatch = formatTaskOutputMatch(task, query);
  if (task.status === "failed" && task.error) {
    return `${task.id} · ${duration} · ${activity}${outputMatch} · failed: ${formatErrorSummary(task.error)}`;
  }
  return `${task.id} · ${duration} · ${activity}${outputMatch} · updated ${updated}`;
}

function formatTaskActivityCounts(task) {
  const messageCount = Number(task.messageCount || 0);
  const traceCount = Number(task.traceCount || 0);
  const resultFileCount = Number(task.resultFileCount || 0);
  return [
    `${messageCount} ${messageCount === 1 ? "msg" : "msgs"}`,
    `${traceCount} ${traceCount === 1 ? "tool" : "tools"}`,
    `${resultFileCount} ${resultFileCount === 1 ? "file" : "files"}`,
  ].join(" · ");
}

function formatSnapshotFailure(snapshot) {
  if (snapshot?.status !== "failed") return "";
  const finalAnswer = snapshot.messages?.findLast((message) => message.role === "assistant")?.content || "";
  return formatErrorSummary(snapshot.error || finalAnswer);
}

function formatErrorSummary(value) {
  const text = String(value || "").trim();
  if (!text) return "";
  const singleQuotedMessage = text.match(/message'\s*:\s*'([^']+)'/);
  if (singleQuotedMessage?.[1]) return singleQuotedMessage[1];
  const doubleQuotedMessage = text.match(/"message"\s*:\s*"([^"]+)"/);
  if (doubleQuotedMessage?.[1]) return doubleQuotedMessage[1];
  const runFailed = text.match(/^Run failed:\s*(.+)$/);
  return runFailed?.[1] || text;
}

async function openTaskFromList(taskId) {
  try {
    await selectTask(taskId);
  } catch (error) {
    restoreRouteAfterNavigationFailure();
    showError(error.message);
  }
}

function taskSearchText(task) {
  return [
    task.title,
    task.task,
    task.id,
    task.status,
    formatStatusLabel(task.status),
    task.error,
    task.dataPath,
    formatTaskMeta(task),
  ].join(" ").toLowerCase();
}

function formatTaskMeta(task) {
  const updated = new Date(task.updatedAt).toLocaleString();
  const dataPath = task.dataPath ? ` · ${task.dataPath}` : "";
  const duration = formatRunDuration(task);
  const activity = formatTaskActivityCounts(task);
  const outputs = formatTaskOutputNames(task);
  if (task.status === "failed" && task.error) {
    return `${task.id}${dataPath} · ${duration} · ${activity}${outputs} · failed: ${formatErrorSummary(task.error)}`;
  }
  return `${task.id}${dataPath} · ${duration} · ${activity}${outputs} · updated ${updated}`;
}

function formatTaskOutputNames(task) {
  const names = task.resultFileNames || [];
  return names.length ? ` · outputs: ${names.join(", ")}` : "";
}

function formatTaskOutputMatch(task, query) {
  const normalizedQuery = String(query || "").trim().toLowerCase();
  if (!normalizedQuery) return "";
  const matches = (task.resultFileNames || [])
    .filter((name) => String(name).toLowerCase().includes(normalizedQuery))
    .slice(0, 3);
  return matches.length ? " · match: " + matches.join(", ") : "";
}

function formatTaskTooltip(task) {
  return `${task.task || task.title || "Untitled BioAgent task"}\n${formatTaskMeta(task)}`;
}

function formatTaskAriaLabel(task) {
  return formatTaskTooltip(task).replaceAll("\n", " · ");
}

function statusClass(status) {
  return `status-${String(status || "idle").replaceAll("_", "-")}`;
}

function formatStatusLabel(value) {
  const status = String(value || "idle");
  if (status === "needs_user_input") {
    return "Needs input";
  }
  return humanizeToolName(status);
}

function isActivationKey(event) {
  return event.key === "Enter" || event.key === " ";
}

function enableKeyboardActivation(row) {
  row.tabIndex = 0;
  row.setAttribute("role", "button");
  row.addEventListener("keydown", (event) => {
    if (event.target !== row) return;
    if (!isActivationKey(event)) return;
    event.preventDefault();
    row.click();
  });
}

function beginSubmit(mode) {
  state.submitToken += 1;
  state.submitMode = mode;
  state.isSubmitting = true;
  renderComposer(state.snapshot);
  return state.submitToken;
}

function cancelSubmit() {
  if (!state.isSubmitting && !state.submitMode) return;
  state.submitToken += 1;
  state.isSubmitting = false;
  state.submitMode = null;
}

function finishSubmit(submitToken) {
  if (state.submitToken !== submitToken) return;
  state.isSubmitting = false;
  state.submitMode = null;
  renderComposer(state.snapshot);
}

function renderComposer(snapshot) {
  if (state.isSubmitting) {
    els.sendButton.textContent = state.submitMode === "reply" ? "Replying" : "Creating";
    els.sendButton.disabled = true;
    els.taskInput.placeholder = state.submitMode === "reply" ? "Sending reply..." : "Creating a new BioAgent run...";
    els.composerOptions.hidden = state.submitMode === "reply";
    return;
  }
  if (snapshot?.status === "needs_user_input") {
    const replyPrompt = latestAssistantPrompt(snapshot);
    els.sendButton.textContent = "Reply";
    els.sendButton.disabled = false;
    els.taskInput.placeholder = replyPrompt ? `Reply to: ${replyPrompt}` : "Reply to BioAgent and continue this run...";
    els.composerOptions.hidden = true;
    if (state.replyTaskId !== snapshot.id) {
      els.taskInput.value = "";
      state.replyTaskId = snapshot.id;
    }
    return;
  }
  state.replyTaskId = null;
  if (isActiveStatus(snapshot?.status)) {
    els.sendButton.textContent = "Running";
    els.sendButton.disabled = true;
    els.taskInput.placeholder = "BioAgent is running. Watch the stream or open a new task.";
    els.composerOptions.hidden = true;
    return;
  }
  if (isTerminalStatus(snapshot?.status)) {
    els.sendButton.textContent = "Rerun";
    els.sendButton.disabled = false;
    els.taskInput.placeholder = "Edit this task or rerun it as a new BioAgent run...";
    els.composerOptions.hidden = false;
    return;
  }
  els.sendButton.textContent = "Run";
  els.sendButton.disabled = false;
  els.taskInput.placeholder = "Describe an analysis task...";
  els.composerOptions.hidden = false;
}

function latestAssistantPrompt(snapshot) {
  const content = snapshot?.messages?.findLast((message) => message.role === "assistant")?.content || "";
  const firstLine = content.trim().split(/\r?\n/)[0].trim();
  return firstLine.length > 96 ? `${firstLine.slice(0, 96)}...` : firstLine;
}

function syncComposerFromSnapshot(snapshot) {
  if (!snapshot || isActiveStatus(snapshot.status) || snapshot.status === "needs_user_input") return;
  els.taskInput.value = snapshot.task || "";
  els.dataPathInput.value = snapshot.dataPath || "";
  els.maxTurnsInput.value = String(snapshot.maxTurns || state.defaultMaxTurns);
}

function renderMessages(snapshot) {
  if (!snapshot) {
    els.messageStream.innerHTML = `<article class="message"><div class="speaker">BioAgent</div><p>Start a task to watch the agent work.</p></article>`;
    setScrollLatestVisible(false);
    return;
  }
  const stickToBottom = shouldStickToBottom(els.messageStream);
  const messages = [{ role: "user", content: snapshot.task }, ...(snapshot.messages || [])];
  els.messageStream.innerHTML = messages.map((message) => `
    <article class="message ${message.role === "user" ? "user" : ""} ${message.error ? "error" : ""} ${message.final ? "final" : ""}">
      <div class="speaker">${message.role === "user" ? "User" : message.final ? "Final" : "BioAgent"}</div>
      <p>${escapeHtml(message.content || "")}</p>
    </article>
  `).join("") + traceBlock(snapshot.traces || [], snapshot);
  if (stickToBottom) {
    scrollToLatest();
  } else {
    setScrollLatestVisible(true);
  }
}

function shouldStickToBottom(element) {
  return element.scrollHeight - element.scrollTop - element.clientHeight < 80;
}

function scrollToLatest() {
  els.messageStream.scrollTop = els.messageStream.scrollHeight;
  setScrollLatestVisible(false);
}

function setScrollLatestVisible(visible) {
  els.scrollLatestButton.classList.toggle("visible", Boolean(visible));
}

function traceBlock(traces, snapshot) {
  if (!traces.length) return "";
  return `
    <section class="trace-block">
      <div class="trace-head"><span>Realtime traces</span><span class="meta">click a row to inspect</span></div>
      ${traces.map((trace) => {
        const displayStatus = formatTraceMeta(trace, snapshot);
        const statusClass = formatTraceStatusClass(trace, snapshot);
        return `
          <div class="trace-row ${trace.id === state.selectedTraceId ? "selected" : ""}" data-trace-id="${trace.id}">
            <span class="status-dot ${escapeHtml(statusClass)}"></span>
            <span>${escapeHtml(formatToolName(trace.toolName))}</span>
            <span class="meta">${escapeHtml(displayStatus)}</span>
          </div>
        `;
      }).join("")}
    </section>
  `;
}

function renderTraces(snapshot) {
  renderTraceBrowser(snapshot);
  const rows = els.messageStream.querySelectorAll("[data-trace-id]");
  rows.forEach((row) => {
    row.addEventListener("click", () => {
      selectTrace(row.dataset.traceId);
    });
    enableKeyboardActivation(row);
  });
  markSelectedTraceRow();
  renderTraceDetail();
}

function renderTraceBrowser(snapshot) {
  const traces = snapshot?.traces || [];
  els.traceListCount.textContent = String(traces.length);
  els.traceList.innerHTML = traces.length ? traces.map((trace, index) => {
    const displayStatus = formatTraceMeta(trace, snapshot);
    const statusClass = formatTraceStatusClass(trace, snapshot);
    return `
      <button class="trace-list-row ${trace.id === state.selectedTraceId ? "selected" : ""}" type="button" data-trace-id="${escapeHtml(trace.id)}" title="${escapeHtml(`${formatToolName(trace.toolName)} · ${displayStatus}`)}">
        <span class="trace-sequence">${String(index + 1).padStart(2, "0")}</span>
        <span class="trace-list-copy">
          <strong>${escapeHtml(formatToolName(trace.toolName))}</strong>
          <span><i class="status-dot ${escapeHtml(statusClass)}"></i>${escapeHtml(displayStatus)}</span>
        </span>
      </button>
    `;
  }).join("") : `<div class="empty">Tool calls appear as the run progresses.</div>`;
  els.traceList.querySelectorAll("[data-trace-id]").forEach((row) => {
    row.addEventListener("click", () => selectTrace(row.dataset.traceId));
  });
  revealSelectedTraceRow();
}

function revealSelectedTraceRow() {
  const selectedRow = [...els.traceList.querySelectorAll("[data-trace-id]")]
    .find((row) => row.dataset.traceId === state.selectedTraceId);
  if (!selectedRow || !els.traceList.clientHeight) return;
  els.traceList.scrollTop = Math.max(0, selectedRow.offsetTop + selectedRow.offsetHeight - els.traceList.clientHeight);
}

function selectTrace(traceId) {
  state.selectedTraceId = traceId;
  markSelectedTraceRow();
  renderTraceDetail();
  activateTab("trace", { user: true });
}

function markSelectedTraceRow() {
  [...els.messageStream.querySelectorAll("[data-trace-id]"), ...els.traceList.querySelectorAll("[data-trace-id]")].forEach((row) => {
    row.classList.toggle("selected", row.dataset.traceId === state.selectedTraceId);
  });
}

function renderTraceDetail() {
  const trace = state.snapshot?.traces?.find((item) => item.id === state.selectedTraceId) || state.snapshot?.traces?.[0];
  if (!trace) {
    els.traceTitle.textContent = "No tool selected";
    els.traceMeta.textContent = "Click a trace row to inspect tool input and output.";
    els.traceInput.textContent = "{}";
    els.traceOutput.textContent = "{}";
    return;
  }
  state.selectedTraceId = trace.id;
  els.traceTitle.textContent = formatToolName(trace.toolName);
  els.traceMeta.textContent = `${formatTraceMeta(trace, state.snapshot)} · ${trace.id}`;
  const missingOutput = trace.status === "running"
    ? "Waiting for tool output..."
    : "No output was received before the call ended.";
  els.traceInput.textContent = trace.input == null
    ? "No input was recorded for this call."
    : formatTracePayload(trace.input, ["code"]);
  els.traceOutput.textContent = trace.output == null
    ? missingOutput
    : formatTracePayload(trace.output, ["stdout", "stderr"]);
}

function formatTraceMeta(trace, snapshot) {
  const status = formatTraceStatus(trace, snapshot);
  const duration = formatRunDuration(trace);
  return `${status} · ${duration}`;
}

function formatTraceStatus(trace, snapshot) {
  if (snapshot?.status === "completed" && trace.status === "failed") {
    return "recovered";
  }
  return formatStatusLabel(trace.status);
}

function formatTraceStatusClass(trace, snapshot) {
  if (snapshot?.status === "completed" && trace.status === "failed") {
    return "recovered";
  }
  return String(trace.status || "pending").replaceAll("_", "-").toLowerCase();
}

function renderPlan(snapshot) {
  const plan = snapshot?.plan || [];
  const visiblePlan = visiblePlanItems(plan);
  els.planCount.textContent = String(visiblePlan.length);
  els.tracePlanCount.textContent = String(visiblePlan.length);
  renderPlanItems(els.planList, visiblePlan, snapshot);
  renderPlanItems(els.tracePlanList, visiblePlan, snapshot);
}

function renderPlanItems(container, visiblePlan, snapshot) {
  container.innerHTML = visiblePlan.length ? "" : `<div class="empty">Plan updates appear as the run progresses.</div>`;
  for (const item of visiblePlan) {
    const step = document.createElement("div");
    const displayStatus = formatPlanStatus(item, snapshot);
    step.className = "step";
    step.innerHTML = `
      <div class="step-line">
        <span class="check ${escapeHtml(displayStatus)}"></span>
        <strong>${escapeHtml(formatPlanTitle(item.title))}</strong>
        <span class="meta">${escapeHtml(displayStatus)}</span>
      </div>
    `;
    container.appendChild(step);
  }
}

function formatPlanStatus(item, snapshot) {
  if (snapshot?.status === "completed" && item.status === "failed") {
    return "recovered";
  }
  return item.status || "pending";
}

function visiblePlanItems(plan) {
  const callTitles = new Set(plan.map((item) => String(item.title || "")).filter((title) => title.startsWith("Call ")));
  return plan.filter((item) => {
    const title = String(item.title || "");
    if (!title.startsWith("Finish ")) return true;
    const toolName = title.slice("Finish ".length);
    return !callTitles.has(`Call ${toolName}`);
  });
}

function formatPlanTitle(value) {
  const title = String(value || "");
  if (title.startsWith("Call ")) {
    return formatToolName(title.slice("Call ".length));
  }
  return title;
}

function formatToolName(value) {
  return humanizeToolName(value);
}

function humanizeToolName(value) {
  const label = String(value || "").replaceAll("_", " ");
  return label.charAt(0).toUpperCase() + label.slice(1);
}

function renderFiles(snapshot) {
  const sourceFiles = snapshot?.resultFiles || [];
  const query = els.fileSearch.value.trim().toLowerCase();
  syncFileSearchClearButton(query);
  const files = filterFileNodes(sourceFiles, query);
  els.fileTree.innerHTML = files.length ? "" : `<div class="empty">${sourceFiles.length ? "No matching files." : "No output files yet."}</div>`;
  for (const node of files) {
    renderFileNode(node, 0);
  }
}

function filterFileNodes(nodes, query) {
  if (!query) return nodes;
  const matches = [];
  for (const node of nodes) {
    const children = filterFileNodes(node.children || [], query);
    if (fileSearchText(node).includes(query) || children.length) {
      matches.push({ ...node, children });
    }
  }
  return matches;
}

function fileSearchText(node) {
  return [
    node.name,
    node.path,
    node.type,
    node.kind,
  ].join(" ").toLowerCase();
}

function renderFileNode(node, depth) {
  const row = document.createElement("div");
  row.className = `file-node ${node.type === "file" ? "is-file" : "is-directory"} depth-${Math.min(depth, 2)} ${node.path === state.selectedFilePath ? "selected" : ""}`;
  if (node.path) {
    row.dataset.filePath = node.path;
    row.title = formatFileDisplayPath(node);
  }
  row.innerHTML = `
    <span>${formatFileKindLabel(node)}</span>
    <strong>${escapeHtml(node.name)}</strong>
    <span>${node.type === "directory" ? (node.children || []).length : formatBytes(node.size || 0)}</span>
  `;
  if (node.type === "file") {
    row.addEventListener("click", () => openFile(node));
    enableKeyboardActivation(row);
  }
  els.fileTree.appendChild(row);
  for (const child of node.children || []) {
    renderFileNode(child, depth + 1);
  }
}

function formatFileKindLabel(node) {
  if (node.type === "directory") return "DIR";
  if (node.kind === "image") return "IMG";
  if (node.kind === "text") return "TXT";
  return "BIN";
}

function formatFileDisplayPath(node) {
  const rawPath = String(node?.path || "");
  const resultRoot = String(state.snapshot?.compute?.resultRoot || "");
  const normalized = rawPath.replaceAll("\\", "/");
  const normalizedRoot = resultRoot.replaceAll("\\", "/").replace(/\/$/, "");
  if (normalizedRoot && normalized.startsWith(`${normalizedRoot}/`)) {
    return normalized.slice(normalizedRoot.length + 1);
  }
  return rawPath;
}

async function openFile(node) {
  activateTab("file", { user: true });
  state.selectedFilePath = node.path;
  markSelectedFileRow();
  els.fileTitle.textContent = node.name;
  els.fileMeta.textContent = formatFileDisplayPath(node);
  els.copyFilePathButton.dataset.copyValue = node.path;
  els.copyFilePathButton.classList.remove("disabled");
  els.copyFilePathButton.setAttribute("aria-disabled", "false");
  els.copyFilePathButton.disabled = false;
  els.downloadLink.href = `/api/tasks/${state.taskId}/files/download?path=${encodeURIComponent(node.path)}`;
  els.downloadLink.classList.remove("disabled");
  els.downloadLink.setAttribute("aria-disabled", "false");
  if (node.kind === "image") {
    els.filePreview.innerHTML = `<img alt="${escapeHtml(node.name)}" src="${els.downloadLink.href}" />`;
    return;
  }
  if (node.kind === "text") {
    const previewTaskId = state.taskId;
    const previewPath = node.path;
    els.filePreview.innerHTML = `<div class="empty">Loading preview...</div>`;
    try {
      const text = await fetch(`/api/tasks/${previewTaskId}/files/content?path=${encodeURIComponent(previewPath)}`).then((r) => {
        if (!r.ok) throw new Error(r.statusText || "Preview request failed");
        return r.text();
      });
      if (state.taskId !== previewTaskId || state.selectedFilePath !== previewPath) return;
      els.filePreview.innerHTML = `<pre>${escapeHtml(text)}</pre>`;
    } catch (error) {
      if (state.taskId !== previewTaskId || state.selectedFilePath !== previewPath) return;
      els.filePreview.innerHTML = `<div class="empty">Unable to load preview: ${escapeHtml(error.message || error)}</div>`;
    }
    return;
  }
  els.filePreview.innerHTML = `<div class="empty">Binary preview is not available. Size: ${formatBytes(node.size || 0)}. Use Download.</div>`;
}

function markSelectedFileRow() {
  document.querySelectorAll(".file-node[data-file-path], .output-file-row[data-file-path]").forEach((row) => {
    row.classList.toggle("selected", row.dataset.filePath === state.selectedFilePath);
  });
}

function resetFilePreview() {
  state.selectedFilePath = null;
  els.fileTitle.textContent = "No file selected";
  els.fileMeta.textContent = "Open an output file to preview text/images or download binary results.";
  els.copyFilePathButton.classList.add("disabled");
  els.copyFilePathButton.setAttribute("aria-disabled", "true");
  els.copyFilePathButton.disabled = true;
  delete els.copyFilePathButton.dataset.copyValue;
  els.downloadLink.href = "#";
  els.downloadLink.classList.add("disabled");
  els.downloadLink.setAttribute("aria-disabled", "true");
  els.filePreview.innerHTML = "";
}

function resetFileSearch() {
  els.fileSearch.value = "";
  syncFileSearchClearButton("");
}

function resetComposerDefaults() {
  els.taskInput.value = state.defaultTaskText;
  els.dataPathInput.value = state.defaultDataPath;
  els.maxTurnsInput.value = state.defaultMaxTurns;
}

function renderResultVisuals(snapshot) {
  const files = flattenFiles(snapshot?.resultFiles || []);
  const images = files.filter((file) => file.kind === "image");
  const summaries = files.filter((file) => file.name.endsWith(".json"));
  const sortedSummaries = summaries.sort(sortResultCards);
  const resultCards = [...images, ...sortedSummaries.slice(0, 2)].sort(sortResultCards);
  els.resultVisuals.innerHTML = "";
  if (snapshot?.status === "failed" && !files.length) {
    return;
  }
  if (!resultCards.length) {
    const emptyMessage = files.length ? "No previewable result cards." : "No result files yet.";
    els.resultVisuals.innerHTML = `<div class="empty">${emptyMessage}</div>`;
    return;
  }
  for (const file of resultCards) {
    if (file.name.endsWith(".json")) {
      const summary = file;
      const summaryTaskId = state.taskId;
      const summaryPath = summary.path;
      const card = document.createElement("div");
      card.className = "summary-card";
      card.title = formatFileDisplayPath(summary);
      card.innerHTML = `<h3>${escapeHtml(summary.name)}</h3><pre>Loading...</pre>`;
      card.addEventListener("click", () => openFile(summary));
      enableKeyboardActivation(card);
      els.resultVisuals.appendChild(card);
      fetch(`/api/tasks/${summaryTaskId}/files/content?path=${encodeURIComponent(summaryPath)}`)
        .then((r) => {
          if (!r.ok) throw new Error(r.statusText || "Preview request failed");
          return r.text();
        })
        .then((text) => {
          if (state.taskId !== summaryTaskId) return;
          if (!flattenFiles(state.snapshot?.resultFiles || []).some((file) => file.path === summaryPath)) return;
          card.querySelector("pre").textContent = formatSummaryCardPreview(text);
        })
        .catch(() => {
          if (state.taskId !== summaryTaskId) return;
          if (!flattenFiles(state.snapshot?.resultFiles || []).some((file) => file.path === summaryPath)) return;
          card.querySelector("pre").textContent = "Unable to load preview.";
        });
      continue;
    }
    if (file.kind === "image") {
      const image = file;
      const card = document.createElement("div");
      card.className = "plot-card";
      card.title = formatFileDisplayPath(image);
      card.innerHTML = `<h3>${escapeHtml(image.name)}</h3><img alt="${escapeHtml(image.name)}" src="/api/tasks/${state.taskId}/files/download?path=${encodeURIComponent(image.path)}" />`;
      card.addEventListener("click", () => openFile(image));
      enableKeyboardActivation(card);
      els.resultVisuals.appendChild(card);
    }
  }
}

function renderResultFileList(snapshot) {
  const files = flattenFiles(snapshot?.resultFiles || []);
  els.resultFileList.innerHTML = files.length ? "" : `<div class="empty">No output files yet.</div>`;
  for (const file of files) {
    const row = document.createElement("div");
    row.className = `output-file-row ${file.path === state.selectedFilePath ? "selected" : ""}`;
    row.dataset.filePath = file.path;
    row.title = `${formatFileDisplayPath(file)} · ${formatBytes(file.size || 0)}`;
    row.innerHTML = `
      <span>${formatFileKindLabel(file)}</span>
      <strong>${escapeHtml(formatFileDisplayPath(file))}</strong>
      <em>${escapeHtml(formatBytes(file.size || 0))}</em>
    `;
    row.addEventListener("click", () => openFile(file));
    enableKeyboardActivation(row);
    els.resultFileList.appendChild(row);
  }
}

function sortResultCards(left, right) {
  return resultCardRank(left) - resultCardRank(right) || left.name.localeCompare(right.name);
}

function resultCardRank(file) {
  const name = String(file.name || "").toLowerCase();
  if (name.includes("umap")) return 0;
  if (name.includes("qc")) return 1;
  if (name.includes("pca")) return 2;
  if (name.includes("cluster")) return 3;
  if (name === "summary.json") return 10;
  if (name.includes("summary") && name.endsWith(".json")) return 11;
  return 20;
}

function formatSummaryCardPreview(text) {
  try {
    const pretty = JSON.stringify(JSON.parse(text), null, 2);
    return limitPreview(pretty, 1200);
  } catch (error) {
    return limitPreview(text, 1200);
  }
}

function limitPreview(value, maxChars) {
  const text = String(value || "");
  if (text.length <= maxChars) return text;
  return `${text.slice(0, maxChars).trimEnd()}\n...`;
}

function renderResultSummary(snapshot) {
  const finalAnswer = (snapshot?.finalAnswer || "").trim();
  const summaryFile = findSummaryFile(snapshot);
  if (snapshot?.status === "failed") {
    const error = formatErrorSummary(snapshot.error || finalAnswer);
    const logDisabled = snapshot.compute?.logPath ? "" : " disabled";
    const traceDisabled = topbarTrace(snapshot) ? "" : " disabled";
    els.resultSummary.innerHTML = `
      <div class="error-summary">
        <strong>Run failed</strong>
        <span>${escapeHtml(error)}</span>
        <div class="result-actions">
          <button class="icon-button result-primary" type="button" data-rerun-current>Rerun</button>
          <button class="icon-button" type="button" data-open-log${logDisabled}>Open log</button>
          <button class="icon-button" type="button" data-open-failed-trace${traceDisabled}>Inspect failed trace</button>
          ${finalAnswer ? '<button class="icon-button" type="button" data-copy="resultFinalAnswerText">Copy error details</button>' : ""}
        </div>
      </div>
      ${finalAnswer ? `<pre id="resultFinalAnswerText">${escapeHtml(finalAnswer)}</pre>` : ""}
    `;
    return;
  }
  if (!finalAnswer && !summaryFile) {
    els.resultSummary.innerHTML = `<div class="empty">Final summary appears here when the run finishes.</div>`;
    return;
  }
  const summaryMarkup = summaryFile ? `
    <div class="kpi-source" id="summaryKpiSource"></div>
    <div class="kpi-grid" id="summaryKpis"></div>
  ` : "";
  const finalAnswerMarkup = finalAnswer ? `<div class="summary-actions"><button class="icon-button" type="button" data-copy="resultFinalAnswerText">Copy summary</button></div><pre id="resultFinalAnswerText">${escapeHtml(finalAnswer)}</pre>` : "";
  els.resultSummary.innerHTML = `${summaryMarkup}${finalAnswerMarkup}`;
  renderSummaryKpis(snapshot);
}

function renderSummaryKpis(snapshot) {
  const summary = findSummaryFile(snapshot);
  const target = document.getElementById("summaryKpis");
  const source = document.getElementById("summaryKpiSource");
  if (!summary || !target || !state.taskId) return;
  if (source) source.textContent = `Metrics from ${summary.name}`;
  const taskId = state.taskId;
  fetch(`/api/tasks/${taskId}/files/content?path=${encodeURIComponent(summary.path)}`)
    .then((response) => {
      if (!response.ok) throw new Error(response.statusText || "Metrics request failed");
      return response.json();
    })
    .then((payload) => {
      if (state.taskId !== taskId) return;
      const qc = payload.qc_stats || {};
      const qcMetrics = payload.qc_metrics || {};
      const medianMito = firstPresent(qc.pct_mito || {}, ["median"]) ?? firstPresent(qcMetrics, ["median_pct_mt", "pct_mito_median"]);
      const cards = [
        ["Cells", firstPresent(payload, ["n_cells", "filtered_cells", "cells_final"])],
        ["Genes", payload.n_genes],
        ["HVGs", firstPresent(payload, ["n_hvg", "n_hvgs", "hvg_count"])],
        ["Clusters", payload.n_clusters],
        ["PCA", firstPresent(payload, ["pca_n_comps", "n_pcs"])],
        ["UMAP", payload.umap_n_comps],
        ["Median genes", firstPresent(qc.n_genes || {}, ["median"]) ?? firstPresent(qcMetrics, ["median_n_genes", "n_genes_median", "n_genes_by_counts_median"])],
        ["Median mito", medianMito == null ? null : `${Number(medianMito).toFixed(2)}%`],
      ].filter(([, value]) => value !== undefined && value !== null && value !== "");
      target.innerHTML = cards.map(([label, value]) => `
        <div class="kpi-card"><span>${escapeHtml(label)}</span><strong>${escapeHtml(formatMetric(value))}</strong></div>
      `).join("");
    })
    .catch(() => {
      if (state.taskId !== taskId) return;
      target.innerHTML = "";
      if (source) source.textContent = "Metrics unavailable";
    });
}

function findSummaryFile(snapshot) {
  const files = flattenFiles(snapshot?.resultFiles || []);
  return files.find((file) => file.name === "summary.json")
    || files.find((file) => file.name.includes("summary") && file.name.endsWith(".json"))
    || null;
}

function firstPresent(object, keys) {
  for (const key of keys) {
    const value = object?.[key];
    if (value !== undefined && value !== null && value !== "") {
      return value;
    }
  }
  return null;
}

function formatDurationMs(ms) {
  const seconds = Math.max(0, Math.round(ms / 1000));
  if (seconds < 60) return `${seconds}s`;
  const minutes = Math.floor(seconds / 60);
  const remainingSeconds = seconds % 60;
  if (minutes < 60) {
    return remainingSeconds ? `${minutes}m ${remainingSeconds}s` : `${minutes}m`;
  }
  const hours = Math.floor(minutes / 60);
  const remainingMinutes = minutes % 60;
  return remainingMinutes ? `${hours}h ${remainingMinutes}m` : `${hours}h`;
}

function formatRunDuration(snapshot) {
  const start = Date.parse(snapshot?.createdAt || "");
  const end = snapshot?.status === "running" ? Date.now() : Date.parse(snapshot?.updatedAt || "");
  if (!Number.isFinite(start) || !Number.isFinite(end) || end < start) {
    return "Pending";
  }
  return formatDurationMs(end - start);
}

function formatTurnProgress(snapshot) {
  const currentTurn = snapshot?.compute?.turn;
  const maxTurns = snapshot?.compute?.maxTurns || snapshot?.maxTurns;
  if (currentTurn === undefined || currentTurn === null || currentTurn === "") {
    if (isActiveStatus(snapshot?.status) || snapshot?.status === "queued") {
      return maxTurns ? `0 / ${maxTurns}` : "Pending";
    }
    if (["completed", "failed"].includes(snapshot?.status)) {
      return "Not recorded";
    }
    return "Pending";
  }
  return maxTurns ? `${currentTurn} / ${maxTurns}` : String(currentTurn);
}

function renderCompute(snapshot) {
  const facts = computeFacts(snapshot);
  renderComputeCards(els.computeGrid, facts);
  renderComputeCards(els.resultComputeGrid, facts);
}

function computeFacts(snapshot) {
  const compute = snapshot?.compute || {};
  return [
    ["Run", compute.runId || "Pending"],
    ["Data Path", snapshot?.dataPath || "None"],
    ["Max Turns", snapshot?.maxTurns || "Pending"],
    ["Turns", formatTurnProgress(snapshot)],
    ["Duration", formatRunDuration(snapshot)],
    ["Model", compute.model || "Pending"],
    ["Base URL", compute.baseUrl || "Pending"],
    ["Current Tool", compute.currentTool || "Idle"],
    ["Run Dir", compute.runDir || "Pending"],
    ["Log Path", compute.logPath || "Pending"],
    ["Result Root", compute.resultRoot || "Pending"],
  ];
}

function renderComputeCards(container, facts) {
  container.innerHTML = facts.map(([label, value]) => `
    <div class="compute-card">
      <div class="compute-label-row">
        <span>${escapeHtml(label)}</span>
        <div class="compute-actions">
          ${isCopyableComputeFact(label, value) ? `<button class="icon-button compute-action" type="button" data-copy data-copy-value="${escapeHtml(value)}">Copy</button>` : ""}
          ${label === "Log Path" && value !== "Pending" ? '<button class="icon-button compute-action" type="button" data-open-log>Open</button>' : ""}
        </div>
      </div>
      <strong title="${escapeHtml(value)}">${escapeHtml(value)}</strong>
    </div>
  `).join("");
}

function isCopyableComputeFact(label, value) {
  if (["Data Path", "Run Dir", "Log Path", "Result Root"].includes(label)) {
    return !["", "None", "Pending"].includes(String(value || ""));
  }
  return false;
}

async function openLog() {
  if (!state.taskId) return;
  const logPath = state.snapshot?.compute?.logPath || "";
  if (!logPath) return;
  const logTaskId = state.taskId;
  activateTab("file", { user: true });
  state.selectedFilePath = "__task_log__";
  markSelectedFileRow();
  els.fileTitle.textContent = "Run log";
  els.fileMeta.textContent = logPath;
  els.copyFilePathButton.dataset.copyValue = logPath;
  els.copyFilePathButton.classList.remove("disabled");
  els.copyFilePathButton.setAttribute("aria-disabled", "false");
  els.copyFilePathButton.disabled = false;
  els.downloadLink.href = `/api/tasks/${state.taskId}/log/download`;
  els.downloadLink.classList.remove("disabled");
  els.downloadLink.setAttribute("aria-disabled", "false");
  els.filePreview.innerHTML = `<div class="empty">Loading log...</div>`;
  try {
    const text = await fetch(`/api/tasks/${logTaskId}/log/content`).then((r) => {
      if (!r.ok) throw new Error(r.statusText || "Log request failed");
      return r.text();
    });
    if (state.taskId !== logTaskId || state.selectedFilePath !== "__task_log__") return;
    els.filePreview.innerHTML = `<pre>${escapeHtml(text)}</pre>`;
  } catch (error) {
    if (state.taskId !== logTaskId || state.selectedFilePath !== "__task_log__") return;
    els.filePreview.innerHTML = `<div class="empty">Unable to load log: ${escapeHtml(error.message || error)}</div>`;
  }
}

function flattenFiles(nodes) {
  return nodes.flatMap((node) => node.type === "directory" ? flattenFiles(node.children || []) : [node]);
}

function activateDefaultTab(snapshot) {
  if (state.userSelectedTab) return;
  if (snapshot?.status === "completed" || snapshot?.status === "failed") {
    activateTab("result");
    return;
  }
  activateTab("trace");
}

function handleDetailTabKeydown(event) {
  if (!event.target.matches(".tab")) return;
  if (!["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
  const tabs = [...els.detailTabs.querySelectorAll(".tab")];
  const currentIndex = tabs.indexOf(event.target);
  let nextIndex;
  if (event.key === "Home") {
    nextIndex = 0;
  } else if (event.key === "End") {
    nextIndex = tabs.length - 1;
  } else if (event.key === "ArrowLeft") {
    nextIndex = (currentIndex - 1 + tabs.length) % tabs.length;
  } else {
    nextIndex = (currentIndex + 1) % tabs.length;
  }
  event.preventDefault();
  const nextTab = tabs[nextIndex];
  activateTab(nextTab.dataset.tab, { user: true });
  nextTab.focus();
}

function activateTab(name, options = {}) {
  if (options.user) {
    state.userSelectedTab = true;
  }
  const activeTab = document.querySelector(".tab.active")?.dataset.tab;
  const tabChanged = activeTab !== name;
  document.querySelectorAll(".tab").forEach((tab) => {
    const isActive = tab.dataset.tab === name;
    tab.classList.toggle("active", isActive);
    tab.setAttribute("aria-selected", String(isActive));
    tab.tabIndex = isActive ? 0 : -1;
  });
  document.querySelectorAll(".detail-pane").forEach((pane) => {
    const isActive = pane.id === `${name}Pane`;
    pane.classList.toggle("active", isActive);
    pane.hidden = !isActive;
  });
  if (name === "trace") {
    revealSelectedTraceRow();
  }
  if (tabChanged) {
    els.detailBody.scrollTop = 0;
  }
}

function addUserMessage(content) {
  const article = document.createElement("article");
  article.className = "message user";
  article.innerHTML = `<div class="speaker">User</div><p>${escapeHtml(content)}</p>`;
  els.messageStream.appendChild(article);
}

async function api(path, options) {
  const response = await fetch(path, options);
  if (!response.ok) {
    let detail = response.statusText;
    try {
      const payload = await response.json();
      detail = formatApiErrorDetail(payload.detail || detail);
    } catch (_) {
      // Keep status text.
    }
    throw new Error(detail);
  }
  return response.json();
}

function formatApiErrorDetail(detail) {
  if (Array.isArray(detail)) {
    return detail.map(formatApiErrorDetail).filter(Boolean).join("; ");
  }
  if (detail && typeof detail === "object") {
    const location = Array.isArray(detail.loc) ? detail.loc.join(".") : "";
    const message = detail.msg || detail.message || JSON.stringify(detail);
    return location ? `${location}: ${message}` : message;
  }
  return String(detail || "Request failed");
}

function showError(message) {
  els.errorLine.textContent = message || "";
}

function showCopyFeedback(button, label) {
  button.dataset.copyOriginalLabel = button.dataset.copyOriginalLabel || button.textContent;
  if (button.copyFeedbackTimer) {
    window.clearTimeout(button.copyFeedbackTimer);
  }
  button.textContent = label;
  button.copyFeedbackTimer = window.setTimeout(() => {
    button.textContent = button.dataset.copyOriginalLabel;
    delete button.dataset.copyOriginalLabel;
    button.copyFeedbackTimer = null;
  }, 1200);
}

function pretty(value) {
  return JSON.stringify(value ?? {}, null, 2);
}

function formatTracePayload(value, expandedKeys = []) {
  if (typeof value === "string") return value;
  if (!value || typeof value !== "object" || Array.isArray(value)) return pretty(value);
  const keys = expandedKeys.filter((key) => Object.hasOwn(value, key));
  if (!keys.length) return pretty(value);
  const metadata = Object.fromEntries(Object.entries(value).filter(([key]) => !keys.includes(key)));
  const sections = Object.keys(metadata).length ? [pretty(metadata)] : [];
  keys.forEach((key) => {
    const text = String(value[key] ?? "");
    sections.push(`${key}:\n${text || "(empty)"}`);
  });
  return sections.join("\n\n");
}

function formatBytes(value) {
  if (value < 1024) return `${value} B`;
  if (value < 1024 * 1024) return `${(value / 1024).toFixed(1)} KB`;
  return `${(value / 1024 / 1024).toFixed(1)} MB`;
}

function formatMetric(value) {
  if (typeof value === "number") {
    return Number.isInteger(value) ? value.toLocaleString() : value.toFixed(2);
  }
  return value;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}
