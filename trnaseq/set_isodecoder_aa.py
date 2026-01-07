import re
import argparse
from collections import defaultdict
import os

def process_trna_data(input_file):
    """
    处理tRNA数据，生成isodecoder和氨基酸水平的统计结果
    
    参数:
        input_file: 输入文件名
    
    返回:
        tuple: (isoacceptor_dict, isodecoder_dict, amino_acid_dict)
    """
    
    # 检查文件是否存在
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"输入文件不存在: {input_file}")
    
    # 读取数据
    data = []
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:  # 跳过空行
                continue
                
            parts = line.split()
            if len(parts) < 2:
                print(f"警告: 第{line_num}行格式不正确，跳过: {line}")
                continue
            
            trna_name = parts[0]
            try:
                reads = int(parts[1])
                data.append((trna_name, reads))
            except ValueError:
                print(f"警告: 第{line_num}行reads数不是整数，跳过: {line}")
                continue
    
    if not data:
        raise ValueError("输入文件中没有有效数据")
    
    # 初始化统计字典
    isoacceptor_dict = defaultdict(int)  # 现有的isoacceptor水平
    isodecoder_dict = defaultdict(int)   # isodecoder水平（按反密码子）
    amino_acid_dict = defaultdict(int)   # 氨基酸水平
    
    # 处理每一行数据
    for trna_name, reads in data:
        # 解析tRNA名称
        # 格式: Homo_sapiens_tRNA-Ala-AGC-1-1
        match = re.match(r'Homo_sapiens_tRNA-([A-Za-z]+)-([A-Z]+)-', trna_name)
        
        if match:
            amino_acid = match.group(1)  # 氨基酸，如 Ala
            anticodon = match.group(2)   # 反密码子，如 AGC
            
            # 现有的isoacceptor水平（就是原始数据）
            isoacceptor_dict[trna_name] += reads
            
            # isodecoder水平：相同反密码子的相加
            isodecoder_key = f"Homo_sapiens_tRNA-{amino_acid}-{anticodon}"
            isodecoder_dict[isodecoder_key] += reads
            
            # 氨基酸水平：相同氨基酸的相加
            amino_acid_key = f"Homo_sapiens_tRNA-{amino_acid}"
            amino_acid_dict[amino_acid_key] += reads
        else:
            print(f"警告: 无法解析tRNA名称格式: {trna_name}")
    
    return isoacceptor_dict, isodecoder_dict, amino_acid_dict

def write_output_files(isoacceptor_dict, isodecoder_dict, amino_acid_dict, 
                      output_isodecoder_file, output_amino_acid_file):
    """
    将统计结果写入文件
    """
    
    # 写入isodecoder结果
    with open(output_isodecoder_file, 'w') as f:
        for key in sorted(isodecoder_dict.keys()):
            f.write(f"{key}\t{isodecoder_dict[key]}\n")
    
    # 写入氨基酸结果
    with open(output_amino_acid_file, 'w') as f:
        for key in sorted(amino_acid_dict.keys()):
            f.write(f"{key}\t{amino_acid_dict[key]}\n")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='处理tRNA数据，生成isodecoder和氨基酸水平的统计结果')
    parser.add_argument('input_file', help='输入文件路径（tRNA数据）')
    parser.add_argument('-o', '--output_dir', default='.', help='输出目录路径（默认为当前目录）')
    parser.add_argument('--isodecoder', default='trna_isodecoder.txt', help='isodecoder输出文件名')
    parser.add_argument('--amino_acid', default='trna_amino_acid.txt', help='氨基酸输出文件名')
    
    args = parser.parse_args()
    
    try:
        # 确保输出目录存在
        os.makedirs(args.output_dir, exist_ok=True)
        
        # 构建输出文件完整路径
        output_isodecoder_file = os.path.join(args.output_dir, args.isodecoder)
        output_amino_acid_file = os.path.join(args.output_dir, args.amino_acid)
        
        print(f"输入文件: {args.input_file}")
        print(f"输出目录: {args.output_dir}")
        print(f"isodecoder输出: {output_isodecoder_file}")
        print(f"氨基酸输出: {output_amino_acid_file}")
        print("-" * 50)
        
        # 处理数据
        isoacceptor_dict, isodecoder_dict, amino_acid_dict = process_trna_data(args.input_file)
        
        # 输出结果到文件
        write_output_files(isoacceptor_dict, isodecoder_dict, amino_acid_dict,
                          output_isodecoder_file, output_amino_acid_file)
        
        # 在控制台显示统计信息
        print(f"\n处理完成！统计信息:")
        print(f"总tRNA种类数: {len(isoacceptor_dict)}")
        print(f"isodecoder种类数: {len(isodecoder_dict)}")
        print(f"氨基酸种类数: {len(amino_acid_dict)}")
        
        # 显示总reads数
        total_reads = sum(isoacceptor_dict.values())
        isodecoder_total = sum(isodecoder_dict.values())
        amino_acid_total = sum(amino_acid_dict.values())
        
        print(f"\n总reads数统计:")
        print(f"原始数据总reads: {total_reads}")
        print(f"isodecoder总reads: {isodecoder_total}")
        print(f"氨基酸总reads: {amino_acid_total}")
        
        # 显示前几个结果作为示例
        print(f"\nisodecoder水平前5个结果:")
        for key in sorted(isodecoder_dict.keys())[:5]:
            print(f"  {key}: {isodecoder_dict[key]}")
            
        print(f"\n氨基酸水平前5个结果:")
        for key in sorted(amino_acid_dict.keys())[:5]:
            print(f"  {key}: {amino_acid_dict[key]}")
            
        print(f"\n结果文件已保存:")
        print(f"  isodecoder: {output_isodecoder_file}")
        print(f"  氨基酸水平: {output_amino_acid_file}")
            
    except FileNotFoundError as e:
        print(f"错误: {e}")
    except ValueError as e:
        print(f"错误: {e}")
    except Exception as e:
        print(f"处理过程中出现错误: {e}")

if __name__ == "__main__":
    main()
