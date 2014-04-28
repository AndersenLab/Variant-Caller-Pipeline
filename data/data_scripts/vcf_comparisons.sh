#!/bin/bash

# Begin by uploading relevant strains:

cd /Volumes/PortusTutus/Raw/FASTQ_clean

# mmp strains

scp  fq/BGI1-RET2-AB1-6b0f1-1.fq.gz fq/BGI1-RET2-AB1-10958-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET2b-AB1-d1442-1.fq.gz fq/BGI3-RET2b-AB1-336eb-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET7-CB4854-cc1ab-1.fq.gz fq/BGI1-RET7-CB4854-9f1f9-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET7-CB4854-f296f-1.fq.gz fq/BGI2-RET7-CB4854-0c00c-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET7-CB4854-12a4d-1.fq.gz fq/BGI2-RET7-CB4854-6eee2-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET7a-CB4854-6e4f4-1.fq.gz fq/BGI3-RET7a-CB4854-b4c61-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET2-ED3017-35456-1.fq.gz fq/BGI1-RET2-ED3017-05974-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET2a-ED3017-1f8fe-1.fq.gz fq/BGI3-RET2a-ED3017-8472b-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET7-ED3040-7228f-1.fq.gz fq/BGI1-RET7-ED3040-1d1fa-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET7-ED3040-befbe-1.fq.gz fq/BGI2-RET7-ED3040-f4344-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET7-ED3040-3cf80-1.fq.gz fq/BGI2-RET7-ED3040-63796-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET7a-ED3040-23d80-1.fq.gz fq/BGI3-RET7a-ED3040-6d1f1-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-ED3049-e0754-1.fq.gz fq/BGI1-RET1-ED3049-e8394-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1b-ED3049-ff88c-1.fq.gz fq/BGI3-RET1b-ED3049-22f9b-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET5-ED3052-5c761-1.fq.gz fq/BGI1-RET5-ED3052-0de79-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET5-ED3052-2a451-1.fq.gz fq/BGI2-RET5-ED3052-1683b-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET5-ED3052-d798f-1.fq.gz fq/BGI2-RET5-ED3052-5010e-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET5b-ED3052-5aee0-1.fq.gz fq/BGI3-RET5b-ED3052-86c51-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET6-GXW1-5c89a-1.fq.gz fq/BGI1-RET6-GXW1-a47bd-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET6-GXW1-86159-1.fq.gz fq/BGI2-RET6-GXW1-f1d13-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET6-GXW1-ec84a-1.fq.gz fq/BGI2-RET6-GXW1-74660-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET6b-GXW1-c2310-1.fq.gz fq/BGI3-RET6b-GXW1-33b2d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET3-JU1088-d7571-1.fq.gz fq/BGI1-RET3-JU1088-01543-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET3b-JU1088-4b749-1.fq.gz fq/BGI3-RET3b-JU1088-e54af-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET4-JU1400-1d78f-1.fq.gz fq/BGI1-RET4-JU1400-b52e7-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET4-JU1400-16dc8-1.fq.gz fq/BGI2-RET4-JU1400-0ce0a-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET4-JU1400-7def4-1.fq.gz fq/BGI2-RET4-JU1400-84c24-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET4a-JU1400-70227-1.fq.gz fq/BGI3-RET4a-JU1400-56504-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET7-JU1652-9c6f7-1.fq.gz fq/BGI1-RET7-JU1652-8f8fe-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET7-JU1652-f0142-1.fq.gz fq/BGI2-RET7-JU1652-ef2d3-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET7-JU1652-d79c3-1.fq.gz fq/BGI2-RET7-JU1652-dca86-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET7b-JU1652-096a4-1.fq.gz fq/BGI3-RET7b-JU1652-93e92-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-JU258-31ebc-1.fq.gz fq/BGI1-RET1-JU258-76968-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-JU258-47483-1.fq.gz fq/BGI3-RET1a-JU258-76ed8-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-JU360-88da4-1.fq.gz fq/BGI1-RET1-JU360-8ce09-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1b-JU360-2b3d8-1.fq.gz fq/BGI3-RET1b-JU360-b44e4-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET4-JU394-5e0da-1.fq.gz fq/BGI1-RET4-JU394-b39b6-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET4-JU394-17512-1.fq.gz fq/BGI2-RET4-JU394-a6a9c-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET4b-JU394-0c024-1.fq.gz fq/BGI3-RET4b-JU394-c4e17-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET5-JU397-0bb39-1.fq.gz fq/BGI1-RET5-JU397-54a80-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET5-JU397-b706f-1.fq.gz fq/BGI2-RET5-JU397-3902e-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET5-JU397-0813b-1.fq.gz fq/BGI2-RET5-JU397-59f4e-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET5b-JU397-5e0d4-1.fq.gz fq/BGI3-RET5b-JU397-c0643-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET4-JU642-a66e8-1.fq.gz fq/BGI1-RET4-JU642-c60a7-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET4-JU642-c6556-1.fq.gz fq/BGI2-RET4-JU642-fd848-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET4b-JU642-f798e-1.fq.gz fq/BGI3-RET4b-JU642-8949e-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET2-JU775-c807d-1.fq.gz fq/BGI1-RET2-JU775-c6ffb-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET2a-JU775-63048-1.fq.gz fq/BGI3-RET2a-JU775-ced96-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET5-KR314-28089-1.fq.gz fq/BGI1-RET5-KR314-3c17d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET5-KR314-26761-1.fq.gz fq/BGI2-RET5-KR314-5011d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET5-KR314-9966f-1.fq.gz fq/BGI2-RET5-KR314-a9768-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET5b-KR314-bc823-1.fq.gz fq/BGI3-RET5b-KR314-69550-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET3-LKC34-767c0-1.fq.gz fq/BGI1-RET3-LKC34-cf838-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET3a-LKC34-b62f9-1.fq.gz fq/BGI3-RET3a-LKC34-1523d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET2-MY1-0547b-1.fq.gz fq/BGI1-RET2-MY1-2c69a-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET2b-MY1-2b685-1.fq.gz fq/BGI3-RET2b-MY1-29fb2-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET2-MY16-42063-1.fq.gz fq/BGI1-RET2-MY16-32ca7-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET2a-MY16-8c14c-1.fq.gz fq/BGI3-RET2a-MY16-b2735-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta

# A few extras:

# C4index issue; NIC4 appeared to be improperly named in BGI3
scp  fq/BGI2-RET4-NIC4-8cb1e-1.fq.gz fq/BGI2-RET4-NIC4-7dc3a-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET4a-NIC4-2e8ba-1.fq.gz fq/BGI3-RET4a-NIC4-7400c-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta

# Header issues:
scp  fq/BGI2-RET6-JU1581-f0534-1.fq.gz fq/BGI2-RET6-JU1581-c4591-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI2-RET6-JU1581-ac5ad-1.fq.gz fq/BGI2-RET6-JU1581-fe817-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta

# Ret1, Ret1a
scp  fq/BGI3-RET1a-AF16-2ea4c-1.fq.gz fq/BGI3-RET1a-AF16-af6fc-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-CB4856_CGC-a40e8-1.fq.gz fq/BGI3-RET1a-CB4856_CGC-1790b-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-DL238-d40d2-1.fq.gz fq/BGI3-RET1a-DL238-9916b-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-HK104-48b14-1.fq.gz fq/BGI3-RET1a-HK104-ff101-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-MY23-d9d91-1.fq.gz fq/BGI3-RET1a-MY23-efaf8-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-N2_CGC-31ee6-1.fq.gz fq/BGI3-RET1a-N2_CGC-94948-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI3-RET1a-QX1430-4038a-1.fq.gz fq/BGI3-RET1a-QX1430-0fa8d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-AF16-f1057-1.fq.gz fq/BGI1-RET1-AF16-e1bb9-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-CB4852-bb482-1.fq.gz fq/BGI1-RET1-CB4852-5af7b-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-CB4856_CGC-1854b-1.fq.gz fq/BGI1-RET1-CB4856_CGC-12b2d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-CB4857_CGC-f2e49-1.fq.gz fq/BGI1-RET1-CB4857_CGC-9b080-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-CX11315-010a1-1.fq.gz fq/BGI1-RET1-CX11315-79b54-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-DL238-47de1-1.fq.gz fq/BGI1-RET1-DL238-6476c-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-HK104-89d1c-1.fq.gz fq/BGI1-RET1-HK104-85e53-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-JU1580-e321c-1.fq.gz fq/BGI1-RET1-JU1580-d9aa4-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-JU778-2383d-1.fq.gz fq/BGI1-RET1-JU778-3962c-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-MY23-362e2-1.fq.gz fq/BGI1-RET1-MY23-0c706-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-N2_CGC-e439a-1.fq.gz fq/BGI1-RET1-N2_CGC-e3536-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-NIC232-9940c-1.fq.gz fq/BGI1-RET1-NIC232-dec8d-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-QG557-60682-1.fq.gz fq/BGI1-RET1-QG557-f8de5-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-QX1430-584e1-1.fq.gz fq/BGI1-RET1-QX1430-f7dab-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
scp  fq/BGI1-RET1-WN2021-afccb-1.fq.gz fq/BGI1-RET1-WN2021-d659a-2.fq.gz dec211@dhunni:/lscr2/andersenlab/dec211/data/fasta
