import os
import base64

def create_html_with_base64_image(image_path, output_html_path):
    """
    Reads a PNG image, converts it to base64, and creates an HTML file with the embedded image.

    Args:
        image_path (str): Path to the input PNG image
        output_html_path (str): Path to the output HTML file
    """
    # Check if the image file exists
    if not os.path.exists(image_path):
        print(f"Error: Image file '{image_path}' does not exist.")
        return

    # Read the image file in binary mode
    with open(image_path, 'rb') as image_file:
        # Encode the image as base64
        base64_encoded = base64.b64encode(image_file.read()).decode('utf-8')

    # Create HTML content with the base64-encoded image
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>UMAP Clustering</title>
</head>
<body>
    <h1>Leiden Clustering UMAP</h1>
    <img src="data:image/png;base64,{base64_encoded}" alt="UMAP Visualization">
</body>
</html>"""

    # Write the HTML content to the output file
    with open(output_html_path, 'w') as html_file:
        html_file.write(html_content)

    print(f"Successfully created HTML file with base64-encoded image: {output_html_path}")