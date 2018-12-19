### create the file of clinical information

library(readxl)
path = "../fNIRs/demog12_18_15/Data book for Hepatic disease study (10001-10055).xlsx"
data_id = read_excel(path, sheet = 1,skip = 2)
data_MH = read_excel(path, sheet = 2)
data_ADL= read_excel(path, sheet = 3)
data_CCI= read_excel(path, sheet = 4)
data_lab= read_excel(path, sheet = 5)
data_MMSE= read_excel(path, sheet = 6,skip = 1)
data_CAM= read_excel(path, sheet = 7)
data_DRS98= read_excel(path, sheet = 9)

delbox_data = data.frame(cam = data_id[,24],dsr = data_id[,25],gender = data_id[,6],age = data_id[,7],YoE = data_MH[,5],meld = data_lab[,21], 
                         cci = as.vector(t(data_CCI[21,3:57])), INR = data_lab[,20], delbox = data_id[,12])
colnames(delbox_data) = c("cam", "DRS-98","gender","age","YoE", "meld","cci", "INR","delbox")

delbox_data$cam = ifelse(delbox_data$cam=="positive",yes = 1,no = 0)
delbox_data$gender = ifelse(delbox_data$gender=="M",yes = 1, no = 0)

delbox_data1 = delbox_data[-c(33,51),] ##subject 33 and 51 didnt complete either SAT or VFT test
delbox_data = delbox_data1[-subj_remove,]

saveRDS(delbox_data, 'demoInfo.rds')
