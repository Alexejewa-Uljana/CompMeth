import numpy as np
from PIL import Image
import os
import sys

sys.path.append(os.path.dirname(__file__))

from main import compress_bmp_svd, decompress_svd


def create_test_images():
    test_images = []
    img1 = np.zeros((200, 200), dtype=np.uint8)
    img1[50:150, 50:150] = 255
    img1[80:120, 80:120] = 128
    Image.fromarray(img1, 'L').save('test1.bmp')
    test_images.append('test1.bmp')
    img2 = np.zeros((150, 150), dtype=np.uint8)
    for i in range(150):
        img2[i, :] = int(i * 255 / 149)
    Image.fromarray(img2, 'L').save('test2.bmp')
    test_images.append('test2.bmp')
    img3 = np.random.randint(0, 256, (180, 180), dtype=np.uint8)
    img3[60:90, 60:120] = 200
    img3[100:130, 40:80] = 100
    Image.fromarray(img3, 'L').save('test3.bmp')
    test_images.append('test3.bmp')
    return test_images


def calculate_quality_metrics(original_path, reconstructed_path):
    original = np.array(Image.open(original_path))
    reconstructed = np.array(Image.open(reconstructed_path))
    mse = np.mean((original - reconstructed) ** 2)
    if mse == 0:
        psnr = float('inf')
    else:
        psnr = 20 * np.log10(255.0 / np.sqrt(mse))
    return mse, psnr


def test_compression_ratios():
    test_images = create_test_images()
    N_values = [2, 4, 8]
    results = []
    for test_img in test_images:
        original_size = os.path.getsize(test_img)
        img = Image.open(test_img)
        original_mode = img.mode
        img.close()

        for N in N_values:
            try:
                base_name = os.path.splitext(test_img)[0]
                compressed_file = f"{base_name}_N{N}.svd"
                reconstructed_file = f"{base_name}_N{N}_reconstructed.bmp"
                k_used = compress_bmp_svd(test_img, compressed_file, N)
                decompress_svd(compressed_file, reconstructed_file, original_mode)
                compressed_size = os.path.getsize(compressed_file)
                actual_ratio = original_size / compressed_size
                mse, psnr = calculate_quality_metrics(test_img, reconstructed_file)
                results.append({
                    'image': test_img,
                    'N': N,
                    'k': k_used,
                    'original_size': original_size,
                    'compressed_size': compressed_size,
                    'actual_ratio': actual_ratio,
                    'mse': mse,
                    'psnr': psnr
                })

            except Exception as e:
                results.append({
                    'image': test_img,
                    'N': N,
                    'error': str(e)
                })
    print_summary_table(results)
    successful = len([r for r in results if 'error' not in r])
    total = len(results)
    assert successful == total


def test_edge_cases():
    small_img = np.ones((20, 20), dtype=np.uint8) * 128
    Image.fromarray(small_img, 'L').save('small_test.bmp')
    edge_cases = [
        ('small_test.bmp', 2),
        ('small_test.bmp', 10),
    ]
    for img_path, N in edge_cases:
        try:
            compressed_file = f"edge_case_N{N}.svd"
            reconstructed_file = f"edge_case_N{N}_reconstructed.bmp"
            compress_bmp_svd(img_path, compressed_file, N)
            decompress_svd(compressed_file, reconstructed_file, 'L')
        except Exception:
            pass


def test_file_format():
    import struct
    test_file = 'test1.bmp'
    compressed_file = 'format_test.svd'
    compress_bmp_svd(test_file, compressed_file, N=4)
    with open(compressed_file, 'rb') as f:
        rows, cols, k = struct.unpack('iii', f.read(12))
        expected_u_size = rows * k * 8
        expected_s_size = k * 8
        expected_v_size = cols * k * 8
        u_data = f.read(expected_u_size)
        s_data = f.read(expected_s_size)
        v_data = f.read(expected_v_size)
        assert len(u_data) == expected_u_size
        assert len(s_data) == expected_s_size
        assert len(v_data) == expected_v_size
        assert len(f.read(1)) == 0


def print_summary_table(results):
    print("Результаты:")
    print("Изображение            |  N  |  k  | Исходный | Сжатый | Коэфф. | MSE    | PSNR")
    print("-" * 85)
    successful_results = [r for r in results if 'error' not in r]
    for res in successful_results:
        img_name = os.path.basename(res['image'])
        print(f"{img_name:20} | {res['N']:3} | {res['k']:3} | {res['original_size']:8} | {res['compressed_size']:6} | {res['actual_ratio']:6.1f} | {res['mse']:6.1f} | {res['psnr']:5.1f}")


def cleanup_test_files():
    test_files = [
        'test1.bmp', 'test2.bmp', 'test3.bmp',
        'small_test.bmp', 'format_test.svd'
    ]
    for base in ['test1', 'test2', 'test3', 'small_test', 'edge_case', 'format_test']:
        for ext in ['.svd', '_reconstructed.bmp']:
            test_files.append(f"{base}{ext}")
        for N in [2, 4, 8, 10]:
            test_files.append(f"{base}_N{N}.svd")
            test_files.append(f"{base}_N{N}_reconstructed.bmp")
    for file in set(test_files):
        if os.path.exists(file):
            try:
                os.remove(file)
            except:
                pass

def main():
    try:
        test_compression_ratios()
        test_edge_cases()
        test_file_format()
        response = input("Удалить тестовые файлы? (y/n): ")
        if response.lower() == 'y':
            cleanup_test_files()
    except Exception as e:
        print(f"Ошибка: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()