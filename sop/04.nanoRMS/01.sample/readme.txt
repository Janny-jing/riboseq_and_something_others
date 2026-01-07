##单个样本分析脚本，使用当前文件夹下的：Pseudou_prediction_singlecondition.R
##single.tsv是个空文件，可以按照格式加入已知的tRNA假鸟苷位点文件
Rscript --vanilla Pseudou_prediction_singlecondition.R [options] -f <epinano_file1> (-s <epinano_file2> -t <epinano_file3>) -p single.tsv

##单样本假鸟苷预测结果文件： CTR1226CTR1_predicted_sites_pU.tsv    CTR1226PIC1_predicted_sites_pU.tsv