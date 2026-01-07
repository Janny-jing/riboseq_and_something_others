#!/usr/bin/env python3
"""
tRNA统计脚本
根据samtools idxstat结果生成isodecoder和氨基酸水平的统计文件
"""

import argparse
import os
import sys

def process_trna_files(input_file, output_isodecoder=None, output_aminoacid=None):
    """
    处理tRNA idxstat结果文件，生成isodecoder和氨基酸水平的统计
    
    参数:
    input_file: 输入文件名
    output_isodecoder: isodecoder输出文件名（可选）
    output_aminoacid: 氨基酸输出文件名（可选）
    """
    
    # 如果未指定输出文件名，自动生成
    if output_isodecoder is None:
        base_name = os.path.splitext(input_file)[0]
        output_isodecoder = f"{base_name}_isodecoder.txt"
    
    if output_aminoacid is None:
        base_name = os.path.splitext(input_file)[0]
        output_aminoacid = f"{base_name}_aminoacid.txt"
    
    # 用于存储isodecoder和氨基酸的统计
    isodecoder_dict = {}
    aminoacid_dict = {}
    
    print(f"正在处理文件: {input_file}")
    
    try:
        # 读取输入文件
        with open(input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:  # 跳过空行
                    continue
                    
                parts = line.split()
                if len(parts) < 3:  # 确保有足够的数据列
                    print(f"警告: 第{line_num}行数据列不足，已跳过")
                    continue
                    
                trna_name = parts[0]
                try:
                    reads_count = int(parts[2])  # 第三列是reads数
                except ValueError:
                    print(f"警告: 第{line_num}行reads数格式错误，已跳过")
                    continue
                
                # 解析tRNA名称
                # 格式: Homo_sapiens_tRNA-Ala-AGC-1-1
                name_parts = trna_name.split('-')
                
                if len(name_parts) >= 4:
                    # 提取isodecoder级别: Homo_sapiens_tRNA-Ala-AGC
                    isodecoder_key = '-'.join(name_parts[:4])
                    
                    # 提取氨基酸级别: Homo_sapiens_tRNA-Ala
                    aminoacid_key = '-'.join(name_parts[:3])
                    
                    # 更新isodecoder统计
                    isodecoder_dict[isodecoder_key] = isodecoder_dict.get(isodecoder_key, 0) + reads_count
                    
                    # 更新氨基酸统计
                    aminoacid_dict[aminoacid_key] = aminoacid_dict.get(aminoacid_key, 0) + reads_count
                else:
                    print(f"警告: 第{line_num}行tRNA名称格式异常: {trna_name}")
        
        # 写入isodecoder结果
        with open(output_isodecoder, 'w') as f:
            for isodecoder, count in sorted(isodecoder_dict.items()):
                f.write(f"{isodecoder}\t{count}\n")
        
        # 写入氨基酸结果
        with open(output_aminoacid, 'w') as f:
            for aminoacid, count in sorted(aminoacid_dict.items()):
                f.write(f"{aminoacid}\t{count}\n")
        
        print(f"\n处理完成！")
        print(f"输入文件: {input_file}")
        print(f"isodecoder水平统计: {output_isodecoder} (共 {len(isodecoder_dict)} 个isodecoder)")
        print(f"氨基酸水平统计: {output_aminoacid} (共 {len(aminoacid_dict)} 个氨基酸类型)")
        print(f"总处理行数: {line_num}")
        
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 处理文件时发生异常 - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='tRNA统计脚本 - 根据samtools idxstat结果生成isodecoder和氨基酸水平的统计文件',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
使用示例:
  %(prog)s trna_results.txt
  %(prog)s trna_results.txt -i isodecoder.txt -a aminoacid.txt
  %(prog)s /path/to/trna_data.txt -i ./output/isodecoder_results.txt
        '''
    )
    
    parser.add_argument('input_file', help='输入的tRNA idxstat结果文件')
    parser.add_argument('-i', '--isodecoder', dest='output_isodecoder', 
                       help='isodecoder输出文件名（可选，不指定则自动生成）')
    parser.add_argument('-a', '--aminoacid', dest='output_aminoacid',
                       help='氨基酸输出文件名（可选，不指定则自动生成）')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input_file):
        print(f"错误: 输入文件 {args.input_file} 不存在")
        sys.exit(1)
    
    # 处理文件
    process_trna_files(args.input_file, args.output_isodecoder, args.output_aminoacid)

if __name__ == "__main__":
    main()