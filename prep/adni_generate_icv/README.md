
# Generating ICV values for ADNI

Intracranial volume (ICV) values are manually generated for ADNI subjects, as ADNI does not seem to provide this information.  This step is more complicated than other steps in this project in that it requires images to be download from ADNI and preprocessed.  **Because of this, the precomputed ICV values are provided in `derivatives/adni_icvs.csv` and do not need to be recomputed to rerun analyses**.  The steps below will document how that table can be created:

## Steps

1. Run `prep/adni_build_table.R` once to generate the table of subjects used in this project.
2. Run `generate_icv_table.R` once to generate a list of ADNI participant IDs (PTIDS.txt) which can be used to search ADNI.
3. Run an image search on ADNI:
    - Go to the ADNI database (ida.loni.usc.edu) > Search > Advanced Image Search (beta).
    - Copy the IDs from PTIDs.txt into the "Subject ID" field.
    - Make sure Image > Modality > MRI and Imaging Protocol > Weighting > T1 are checked,
    as well as Image Types > Original and Image Types > Pre-processed (on the side)
    - Add ALL the resulting images to an image collection.
    - Under your image collections, download the CSV record for the collection just created.  Save this record as MRISEARCH.csv.
    - **NOTE: You do not have to download all the images found in this search.**
4. Rerun `generate_icv_table.R`, which should now spit out DOWNLOAD_IDS.txt.
5. Download images from ADNI:
    - Go to the ADNI database (ida.loni.usc.edu) > Search > Advanced Image Search (beta).
    - Copy the contents of DOWNLOAD_IDS.txt into the "Image ID" field
    - Add the search results to a collection and download them in NIFTI format.
    - In this directory, run `find . -name "ADNI*.nii.gz" > DOWNLOADLIST.txt` and copy the output to this folder.
6. Preprocess the downloaded images to generate an ICV mask.  This step requires ANTS software (http://stnava.github.io/ANTs/) and CBICA DeepMRSeg (https://github.com/CBICA/DeepMRSeg); see the respective websites for help with installation.  You can use `process_image.sh`, or use the sub scripts `preskullstripping.sh` and `skullstripping_deepmrseg.sh` and stitch them into your own pipeline.  You can also try to use the `icv_pipeline.sh` script for applying these tools, though **note that this assumes a SLURM compute cluster is being used**.   If using `icv_pipeline.sh`, make sure you update variables:
     - `INPUT_DIRECTORY` should point the the folder where you downloaded all the ADNI images in DOWNLOAD_IDS.txt
     - `TTP2` should point to the folder containing `preskullstripping.sh` and `skullstripping_deepmrseg.sh` (i.e., this folder.  TTP2=Tom T1 Pipeline 2, where the scripts are pulled from).
     - `SCRIPT` should point to `process_image.sh`.
     - `TTP2` and `SCRIPT` may need to be updated to be absolute paths.
7. Once the images are processed and brain masks have been created, use `mean_intensity_tool.py` to generate mask sizes for all the brain masks computed by DeepMRSeg.  Save the output as MASKSIZES.csv.
8. Run `generate_icv_table.R` for a third time, which should now run to completion and generate `derivatives/adni_icvs.csv`.