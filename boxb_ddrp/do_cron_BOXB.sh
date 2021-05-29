#!/bin/sh
# run cron job to create BOXB  maps for CAPS with DDRP v2 for plant diseases
cd /usr/local/dds/DDRP_B1
./v25dis_risk_DDRP_B2.R --spp BOXB --forecast_data PRISM --start_year 2020 --start_doy 1 --end_doy 365 --keep_leap 0 --region_param EAST --exclusions_stressunits 1 --ref_date 20200615 --pems 1 --mapA 1 --mapE 0 --mapL 0 --mapP 0 --out_option 1 --out_dir BOXB_test
