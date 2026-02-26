## COW

python3 PSMC-FAC.py \
  --root 1.ALL_psmc_files_to_tsv/cattle \
  --base_files 801-Cattle.psmc 911-Cattle.psmc 927-Cattle.psmc \
  --outdir test_output/cattle \
  --mu 0.98e-8 \
  --g 5 \
  --FNR_min 0.0 \
  --FNR_max 0.99 \
  --svalue 100 \
  --tmin 5e4 \
  --tmax 1.5e6 \
  --write_psmc_tsv

## HUMAN

python3  PSMC-FAC.py \
  --root 1.ALL_psmc_files_to_tsv/human \
  --base_files 143-Human.psmc 144-Human.psmc 159-Human.psmc 167-Human.psmc 124-Human.psmc 196-Human.psmc 254-Human.psmc 259-Human.psmc 266-Human.psmc \
  --outdir test_output/human \
  --mu 1.25e-8 \
  --g 25 \
  --FNR_min 0.0 \
  --FNR_max 0.99 \
  --svalue 100 \
  --tmin 5e4 \
  --tmax 1.5e6 \
  --write_psmc_tsv

  ## WOLF

  python3 PSMC-FAC.py \
  --root 1.ALL_psmc_files_to_tsv/wolf \
  --base_files gr_wolf.spain.CluJAL7609.psmc gr_wolf.spain.JAL7487.psmc gr_wolf.spain.SAMN04851099.psmc gr_wolf.italy.SAMEA116045429.psmc gr_wolf.italy.SAMEA116045431.psmc gr_wolf.italy.SAMEA116045435.psmc \
  --outdir test_output/wolf \
  --mu 4.5e-9 \
  --g 4.4 \
  --FNR_min 0.0 \
  --FNR_max 0.99 \
  --svalue 100 \
  --tmin 5e4 \
  --tmax 1.5e6 \
  --write_psmc_tsv