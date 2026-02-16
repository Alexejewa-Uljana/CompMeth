import numpy as np
from scipy.linalg import svd
from PIL import Image
import struct
import os


class SVCFile:
    def __init__(self):
        self.magic = b'SVDC'
        self.version = 1
        self.mode = 'L'
        self.rows = 0
        self.cols = 0
        self.k = 0
        self.U = None
        self.S = None
        self.V = None


def create_svc_from_image(img_path, N):
    img = Image.open(img_path)
    if img.mode not in ['1', 'L']:
        raise ValueError("Изображение должно быть монохромным")

    img_array = np.array(img, dtype=np.float32)
    rows, cols = img_array.shape

    if img.mode == '1':
        size_bmp = (rows * cols + 7) // 8
    else:
        size_bmp = rows * cols

    header_size = 24
    k_max = 0

    for k in range(1, min(rows, cols)):
        size_compressed = header_size + (rows * k + k + cols * k) * 4
        if size_compressed <= size_bmp / N:
            k_max = k
        else:
            break

    if k_max == 0:
        raise ValueError(f"N={N} слишком велико, сжатие невозможно")

    U, S, Vt = svd(img_array, full_matrices=False)

    svc = SVCFile()
    svc.mode = img.mode
    svc.rows = rows
    svc.cols = cols
    svc.k = k_max
    svc.U = U[:, :k_max].astype(np.float32)
    svc.S = S[:k_max].astype(np.float32)
    svc.V = Vt[:k_max, :].T.astype(np.float32)

    return svc


def write_svc_file(svc, output_path):
    with open(output_path, 'wb') as f:
        f.write(svc.magic)
        f.write(struct.pack('H', svc.version))
        mode_flag = 0 if svc.mode == '1' else 1
        f.write(struct.pack('H', mode_flag))
        f.write(struct.pack('III', svc.rows, svc.cols, svc.k))
        f.write(svc.U.tobytes())
        f.write(svc.S.tobytes())
        f.write(svc.V.tobytes())


def read_svc_file(input_path):
    with open(input_path, 'rb') as f:
        magic = f.read(4)
        if magic != b'SVDC':
            raise ValueError("Неверный формат файла")

        svc = SVCFile()
        svc.version = struct.unpack('H', f.read(2))[0]
        mode_flag = struct.unpack('H', f.read(2))[0]
        svc.mode = '1' if mode_flag == 0 else 'L'
        svc.rows, svc.cols, svc.k = struct.unpack('III', f.read(12))

        svc.U = np.frombuffer(f.read(svc.rows * svc.k * 4), dtype=np.float32).reshape(svc.rows, svc.k)
        svc.S = np.frombuffer(f.read(svc.k * 4), dtype=np.float32)
        svc.V = np.frombuffer(f.read(svc.cols * svc.k * 4), dtype=np.float32).reshape(svc.cols, svc.k)

    return svc


def compress_bmp_svd(input_path, output_path, N):
    svc = create_svc_from_image(input_path, N)
    write_svc_file(svc, output_path)

    img = Image.open(input_path)
    if img.mode == '1':
        size_bmp = (svc.rows * svc.cols + 7) // 8
    else:
        size_bmp = svc.rows * svc.cols

    final_size = os.path.getsize(output_path)
    compression_ratio = size_bmp / final_size

    print(f"Сжатие: k={svc.k}, сжатие={compression_ratio:.2f}x")

    return svc.k, final_size


def decompress_svd(input_path, output_path):
    svc = read_svc_file(input_path)

    reconstructed = svc.U @ np.diag(svc.S) @ svc.V.T

    if svc.mode == '1':
        reconstructed = np.where(reconstructed > 127, 255, 0).astype(np.uint8)
    else:
        reconstructed = np.clip(reconstructed, 0, 255).astype(np.uint8)

    Image.fromarray(reconstructed, mode=svc.mode).save(output_path)


def create_example_images():
    if not os.path.exists('images'):
        os.makedirs('images')

    img1 = np.zeros((100, 150), dtype=np.uint8)
    img1[30:70, 50:100] = 255
    Image.fromarray(img1, 'L').save('images/example1.bmp')

    img2 = np.zeros((120, 120), dtype=np.uint8)
    for i in range(120):
        img2[i, :] = int(i * 255 / 119)
    Image.fromarray(img2, 'L').save('images/example2.bmp')

    img3 = np.random.randint(0, 255, (80, 80), dtype=np.uint8)
    Image.fromarray(img3, 'L').save('images/example3.bmp')

    return ['images/example1.bmp', 'images/example2.bmp', 'images/example3.bmp']


def process_examples():
    if not os.path.exists('images'):
        os.makedirs('images')

    examples = ['images/example1.bmp', 'images/example2.bmp', 'images/example3.bmp']
    existing_examples = [ex for ex in examples if os.path.exists(ex)]

    if len(existing_examples) < len(examples):
        examples = create_example_images()

    N_values = [2, 4, 8]
    results_summary = []

    for example in examples:
        if not os.path.exists(example):
            continue

        print(f"Обработка: {os.path.basename(example)}")

        for N in N_values:
            print(f"N = {N}")

            base_name = os.path.splitext(os.path.basename(example))[0]
            compressed_file = f"images/{base_name}_N{N}.svc"
            reconstructed_file = f"images/{base_name}_N{N}_reconstructed.bmp"

            try:
                k_used, compressed_size = compress_bmp_svd(example, compressed_file, N)
                decompress_svd(compressed_file, reconstructed_file)

                original_size = os.path.getsize(example)
                ratio = original_size / compressed_size

                original_img = np.array(Image.open(example))
                reconstructed_img = np.array(Image.open(reconstructed_file))
                mse = np.mean((original_img - reconstructed_img) ** 2)
                psnr = 20 * np.log10(255.0 / np.sqrt(mse)) if mse > 0 else float('inf')

                results_summary.append({
                    'image': example,
                    'N': N,
                    'k': k_used,
                    'original_size': original_size,
                    'compressed_size': compressed_size,
                    'ratio': ratio,
                    'mse': mse,
                    'psnr': psnr
                })

            except Exception as e:
                print(f"Ошибка: {e}")
                results_summary.append({
                    'image': example,
                    'N': N,
                    'error': str(e)
                })

    print("\nРезультаты:")
    print("Изображение    |  N  |  k  | Сжатие | MSE  | PSNR")
    print("-" * 50)

    for result in results_summary:
        if 'error' in result:
            print(
                f"{os.path.basename(result['image']):14} | {result['N']:3} | {'ERROR':3} | {'-':6} | {'-':4} | {'-':4}")
        else:
            img_name = os.path.basename(result['image'])
            print(
                f"{img_name:14} | {result['N']:3} | {result['k']:3} | {result['ratio']:6.1f}x | {result['mse']:4.0f} | {result['psnr']:4.0f}")

    with open('images/summary.txt', 'w', encoding='utf-8') as f:
        f.write("Сводка результатов\n")
        for result in results_summary:
            if 'error' not in result:
                img_name = os.path.basename(result['image'])
                f.write(f"{img_name} N={result['N']}: k={result['k']}, сжатие={result['ratio']:.1f}x\n")

    print("Все файлы в 'images/'")


if __name__ == "__main__":
    process_examples()