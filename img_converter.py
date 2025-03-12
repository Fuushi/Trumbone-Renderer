from PIL import Image

def ppm_to_jpeg(ppm_path, jpeg_path):
    # Open the PPM file
    with Image.open(ppm_path) as img:
        # Convert and save as JPEG
        img.convert("RGB").save(jpeg_path, "JPEG")
        print(f"Converted {ppm_path} to {jpeg_path}")

# Example usage
for i in range(60):
    print(i)

    path = f"./tmp/f{i}.ppmcomposite.ppm"
    opath = path.replace(".ppm", ".jpeg")

    ppm_to_jpeg(path, opath)

