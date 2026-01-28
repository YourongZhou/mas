# Notebooks 使用说明

## test_single_cell_analysis.ipynb

### 功能说明

这个 notebook 用于测试完整的单细胞数据分析与细胞标注流程。

### 预期工作流

1. **Supervisor** → 选择 `code_dev` agent
2. **Code Dev Agent** → 生成单细胞分析代码（UMAP、Leiden 聚类、DE 分析）
3. **Critic** → 审核代码
4. **Supervisor** → 选择 `tool_caller` agent
5. **Tool Caller Agent** → 使用细胞类型注释工具，根据 DE 基因确定细胞类型
6. **Critic** → 审核工具调用结果
7. **Supervisor** → 再次选择 `code_dev` agent
8. **Code Dev Agent** → 使用确定的细胞类型进行标注
9. **Finalize** → 生成最终答案

### 使用步骤

1. **配置环境变量**
   ```python
   # 在 notebook 的第一个 cell 中设置
   os.environ["OPENAI_API_KEY"] = "your_api_key"
   ```

2. **准备数据**
   - 确保单细胞数据文件存在（如 `pbmc3k.h5ad`）
   - 更新 notebook 中的数据路径

3. **运行 notebook**
   - 按顺序执行所有 cell
   - 观察执行流程和结果

4. **查看结果**
   - 执行步骤记录在 `execution_steps` 中
   - 生成的代码保存在 `generated_codes` 中
   - 细胞类型注释结果在 `celltype_results` 中
   - 最终结果保存到 `../results/` 目录

### 输出文件

- `results/single_cell_analysis_result.json`: 执行摘要
- `results/generated_code.py`: 所有生成的代码

### 注意事项

- 需要配置 `.env` 文件或设置环境变量
- 需要网络连接以调用 LLM API
- 执行时间取决于 LLM 响应速度和网络状况
- 如果某个步骤失败，可以手动调整状态后继续执行

