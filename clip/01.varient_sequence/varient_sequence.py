from Bio import SeqIO
from Bio.Seq import Seq
import re
from typing import List, Tuple, Dict

class CovidRegions:
    """新冠病毒基因组区域定义"""
    def __init__(self):
        # 定义各个区域的起始和终止位置（1-based）
        self.regions = {
            "5'UTR": (1, 265),
            "ORF1a": (266, 13468),
            "ORF1b": (13468, 21555),
            "S": (21563, 25384),
            "ORF3a": (25393, 26220),
            "E": (26245, 26472),
            "M": (26523, 27191),
            "ORF6": (27202, 27387),
            "ORF7a": (27394, 27759),
            "ORF7b": (27756, 27887),
            "ORF8": (27894, 28259),
            "N": (28274, 29533),
            "ORF10": (29558, 29674),
            "3'UTR": (29675, 29903)
        }
        
        # 定义核糖体移码位点
        self.frameshift_site = 13468

        # 定义NSP蛋白的位置（考虑移码）
        self.nsp_regions = {
            "NSP1": (266, 805),
            "NSP2": (806, 2719),
            "NSP3": (2720, 8554),
            "NSP4": (8555, 10054),
            "NSP5": (10055, 10972),
            "NSP6": (10973, 11842),
            "NSP7": (11843, 12091),
            "NSP8": (12092, 12685),
            "NSP9": (12686, 13024),
            "NSP10": (13025, 13441),
            "NSP11": (13442, 13480),
            "NSP12": (13468, 16236),
            "NSP13": (16237, 18039),
            "NSP14": (18040, 19620),
            "NSP15": (19621, 20658),
            "NSP16": (20659, 21552)
        }
        
        # 合并所有区域
        self.all_regions = {**self.regions, **self.nsp_regions}

    def get_region_bounds(self, region_name: str) -> Tuple[int, int]:
        """获取指定区域的起始和终止位置"""
        region = self.all_regions.get(region_name)
        if not region:
            raise ValueError(f"未找到区域：{region_name}")
        return region

    def list_all_regions(self) -> List[str]:
        """返回所有可用区域的列表"""
        return sorted(self.all_regions.keys())

    def get_region_by_position(self, position: int) -> str:
        """根据位置获取所在区域名称"""
        for name, (start, end) in self.all_regions.items():
            if start <= position <= end:
                return name
        return "区域外"

class CovidMutationGenerator:
    def __init__(self, reference_path):
        # 读取参考序列
        try:
            self.reference = str(SeqIO.read(reference_path, "fasta").seq)
        except Exception as e:
            raise Exception(f"读取参考序列失败：{str(e)}")
        self.regions = CovidRegions()
        
        # 定义密码子表
        self.genetic_code = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }

    def validate_mutation(self, mutation: str) -> Tuple[int, str]:
        """验证突变格式是否正确"""
        pattern = r'^[ATCG]\d+[ATCG]$'
        if not re.match(pattern, mutation):
            raise ValueError("突变格式不正确，应为如'C241T'的格式")
        
        orig_base = mutation[0]
        pos = int(mutation[1:-1]) - 1  # 转换为0-based索引
        new_base = mutation[-1]
        
        if pos >= len(self.reference):
            raise ValueError("位置超出序列长度范围")
        if self.reference[pos] != orig_base:
            raise ValueError(f"参考序列在位置{pos+1}的碱基为{self.reference[pos]}，而不是{orig_base}")
            
        return pos, new_base

    def check_amino_acid_change(self, mutation: str) -> Tuple[bool, str]:
        """检查是否导致氨基酸改变，考虑核糖体移码"""
        pos, new_base = self.validate_mutation(mutation)
        
        # 获取所在区域
        mutation_region = self.regions.get_region_by_position(pos + 1)
        
        # 确定是否在编码区域
        if mutation_region in ["5'UTR", "3'UTR"]:
            return False, f"位于非编码区域：{mutation_region}"
        
        # 确定读框
        if pos >= self.regions.frameshift_site - 1 and mutation_region in ["ORF1b", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16"]:
            # 考虑-1移码后的读框
            frame_offset = -1
            relative_pos = pos - self.regions.frameshift_site + 1
            codon_start = ((relative_pos // 3) * 3) + self.regions.frameshift_site - 1
        else:
            # 正常的读框
            frame_offset = 0
            codon_start = (pos // 3) * 3
        
        # 获取密码子
        orig_codon = self.reference[codon_start:codon_start+3]
        temp_seq = list(self.reference)
        temp_seq[pos] = new_base
        new_codon = ''.join(temp_seq[codon_start:codon_start+3])
        
        # 翻译
        orig_aa = self.genetic_code.get(orig_codon, 'X')
        new_aa = self.genetic_code.get(new_codon, 'X')
        
        # 准备输出信息
        frame_info = f"读框偏移：{frame_offset}" if frame_offset != 0 else "标准读框"
        position_info = f"密码子位置：{codon_start+1}-{codon_start+3}"
        
        if orig_aa == new_aa:
            return False, f"区域：{mutation_region}\n{frame_info}\n{position_info}\n同义突变：密码子从{orig_codon}变为{new_codon}，氨基酸保持为{orig_aa}"
        else:
            return True, f"区域：{mutation_region}\n{frame_info}\n{position_info}\n非同义突变：密码子从{orig_codon}变为{new_codon}，氨基酸从{orig_aa}变为{new_aa}"

    def generate_mutant(self, mutations: List[str], region: str = None) -> str:
        """生成突变序列，考虑核糖体移码"""
        if region:
            start, end = self.regions.get_region_bounds(region)
            # 转换为0-based索引
            start -= 1
            sequence = list(self.reference[start:end])
            
            # 调整突变位置到区域内的相对位置
            for mutation in mutations:
                pos, new_base = self.validate_mutation(mutation)
                if start <= pos < end:
                    if pos >= self.regions.frameshift_site - 1 and region in ["ORF1b", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16"]:
                        print(f"注意：突变 {mutation} 位于核糖体移码位点之后，使用移码后的读框")
                    sequence[pos - start] = new_base
                else:
                    print(f"警告：突变 {mutation} 在所选区域 {region} 之外，已忽略")
        else:
            sequence = list(self.reference)
            for mutation in mutations:
                pos, new_base = self.validate_mutation(mutation)
                sequence[pos] = new_base
                if pos >= self.regions.frameshift_site - 1:
                    mutation_region = self.regions.get_region_by_position(pos + 1)
                    if mutation_region in ["ORF1b", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16"]:
                        print(f"注意：突变 {mutation} 位于核糖体移码位点之后，使用移码后的读框")
                
        return ''.join(sequence)

def save_sequence(sequence: str, filename: str, region: str = None):
    """保存序列到FASTA文件"""
    with open(filename, 'w') as f:
        header = ">mutated_sequence"
        if region:
            header += f"_{region}"
        f.write(header + "\n")
        # 每行写入60个碱基
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\n')

def main():
    # 使用示例
    ref_path = "wuhan_ref.fasta"  # 需要替换为实际的参考序列文件路径
    try:
        generator = CovidMutationGenerator(ref_path)
    except Exception as e:
        print(f"初始化失败：{str(e)}")
        return
    
    mutations = []
    selected_region = None
    
    while True:
        print("\n当前状态：")
        print(f"选定区域：{selected_region if selected_region else '全基因组'}")
        print("当前的突变列表：", ', '.join(mutations) if mutations else "暂无")
        
        print("\n请选择操作：")
        print("1. 选择区域")
        print("2. 添加新的突变")
        print("3. 删除已有突变")
        print("4. 清空突变列表")
        print("5. 生成突变序列")
        print("6. 退出程序")
        
        choice = input("请输入选项（1-6）：")
        
        if choice == '1':
            print("\n可用的区域：")
            regions = generator.regions.list_all_regions()
            for i, region in enumerate(regions):
                print(f"{i+1}. {region}")
            print(f"{len(regions)+1}. 全基因组")
            
            region_choice = input("请选择区域编号（输入0返回）：")
            try:
                region_idx = int(region_choice)
                if region_idx == 0:
                    continue
                elif region_idx == len(regions)+1:
                    selected_region = None
                    print("已选择全基因组")
                elif 1 <= region_idx <= len(regions):
                    selected_region = regions[region_idx-1]
                    print(f"已选择区域：{selected_region}")
                else:
                    print("无效的选择")
            except ValueError:
                print("请输入有效的数字")
                
        elif choice == '2':
            mutation = input("请输入突变（格式如C241T）：")
            try:
                # 验证突变格式
                generator.validate_mutation(mutation)
                # 检查氨基酸变化
                has_aa_change, aa_info = generator.check_amino_acid_change(mutation)
                
                print(f"\n突变信息：{mutation}")
                print(aa_info)
                print(f"是否导致氨基酸改变：{'是' if has_aa_change else '否'}")
                
                # 如果选择了特定区域，检查突变是否在区域内
                if selected_region:
                    start, end = generator.regions.get_region_bounds(selected_region)
                    pos = int(mutation[1:-1])
                    if not (start <= pos <= end):
                        print(f"警告：此突变位点在所选区域 {selected_region} ({start}-{end})之外")
                
                confirm = input("是否添加这个突变？(y/n): ")
                if confirm.lower() == 'y':
                    mutations.append(mutation)
                    print("突变已添加")
                
            except ValueError as e:
                print(f"错误：{str(e)}")
                
        elif choice == '3':
            if not mutations:
                print("没有可删除的突变")
                continue
                
            print("当前突变列表：")
            for i, mut in enumerate(mutations):
                print(f"{i+1}. {mut}")
            idx = input("请输入要删除的突变编号：")
            try:
                idx = int(idx) - 1
                if 0 <= idx < len(mutations):
                    removed = mutations.pop(idx)
                    print(f"已删除突变：{removed}")
                else:
                    print("无效的编号")
            except ValueError:
                print("请输入有效的数字")
                
        elif choice == '4':
            mutations.clear()
            print("已清空突变列表")
            
        elif choice == '5':
            if not mutations:
                print("请先添加突变")
                continue
                
            print("\n当前状态：")
            print(f"选定区域：{selected_region if selected_region else '全基因组'}")
            print("突变列表：", ', '.join(mutations))
            confirm = input("确定要生成突变序列吗？(y/n): ")
            
            if confirm.lower() == 'y':
                try:
                    mutant_seq = generator.generate_mutant(mutations, selected_region)
                    output_file = f"mutated_sequence{'_'+selected_region if selected_region else ''}.fasta"
                    save_sequence(mutant_seq, output_file, selected_region)
                    print(f"\n突变序列已生成并保存到文件：{output_file}")
                except Exception as e:
                    print(f"生成序列时发生错误：{str(e)}")
            
        elif choice == '6':
            break
            
        else:
            print("无效的选项，请重新选择")

if __name__ == "__main__":
    main()