# 2026-stickleback-intraguild-predation
Data and code from Fassora, J., Martin, J. S., Matthews, B. "Multivariate morphological divergence due to intraguild predation", _Evolution_, **2026**. Study on morphological divergence in Greenlandic sticklebacks due to intraguild predation by Arctic charr.

**Contents of this repository**
- stickle_igp_morphology.Rproj: file to open the R project
- src/trait_analysis_pub.r: script containing the full analysis and generating the plots represented in the paper
- src/mvtraitmod1.stan: STAN model for the covariance matrix that ignores charr presence (called by trait_analysis_pub.r)
- src/smvtraitmod2.stan: STAN model for the covariance matrix that accounts for charr presence (called by trait_analysis_pub.r)
- data/fish_data_pub.csv: data of all measured morphological traits for each fish (for the raw landmarks coordinates contact the authors)
- data/lake_info_pub.csv: data on each sampled lake
- model_outputs: folder that contains a split zip archive with the exact model outputs analysed and discussed in the paper. Note that all parts of the archive are required to extract the files

**Data description**

Fish data (fish_data_pub.csv)
- FishEc: unique fish ID
- Site: lake where the fish was sampled (internal name, see trait_analysis_pub.r for the "translation" to the names displayed in the figures)
- Year: sampling year (note that all fish from a given site were sampled in the same year)
- SLF: "Standard Length Fish", the length from the tip of the mouth to the caudal peduncule (i.e. the tail without tail fins)
- Other columns: measure in millimetres of each trait, see Figure 2 along with Tables S2 and S3 (in the supplementary material) to see what each measure exactly is

Lake data (lake_info_pub.csv)
- Site: unique lake name (internal name, since this can appear confusing and without logic, for the paper and the figures we used an alternative nomenclature, the translation is in the first few lines of trait_analysis_pub.r)
- Char: presence (1) or absence (0) of Arctic charr
- Drainage: factor for the drainage basin of the lake (the names of the drainages have no particular meaning)
- Region: factor for the geographical region of the lake (A: Akia island, T: Tuttutooq island, Q: Qassiarsuk region)
- Lat: latitude in deg
- Lon: longitude in deg
- Area_ha: area of the lake in hectares
- Perimeter_m: perimeter of the lake in metres
