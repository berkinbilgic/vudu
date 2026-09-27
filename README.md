# vudu
## script_vudu_gre.m: 
Performs reconstruction for gradient echo EPI data acquired with Fleet slice ordering using variable flip angles. 
Shot 1 is acquired with AP, and shot 2 is acquired with PA phase encoding.
Separate Sense reconstruction for each shot is compared against Hankel low-rank regularized joint reconstruction in VUDU.

## Data download

The example data (4 files, 429 MB; three are larger than GitHub's 100 MB file limit) are attached to the [v1.0 release](https://github.com/berkinbilgic/vudu/releases/tag/v1.0). Download them into a `data/` folder in the repository:

```bash
mkdir -p data && cd data
for f in Img_epi_ap.mat Img_epi_pa.mat receive.mat img_fieldmap.mat; do curl -LO https://github.com/berkinbilgic/vudu/releases/download/v1.0/$f; done
```

Then run `script_vudu_gre.m` from the repository folder in MATLAB.
