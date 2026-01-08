# ProtruJ_Morphometrics

* **Developed by:** Thomas 
* **Developed for:** Laurent & Emmanuelle
* **Team:** Muller
* **Date:** January 2026
* **Software:** Fiji


### Images description

3D images taken with x10 objective on a spinning-disk microscope.

3 channels:
  1. DAPI nuclei (not used)
  2. GFP protrusions
  3. Brightfield pores

### Macro description

* Detect pores in the 3rd channel using Multi-Template-Matching
* Segment protrusions in MIP 2nd channel 
* For each detected pore:
      * Use the 3D watershed Split algorithm and detect the anastomosis between pores
      * Run skeletonization filter branches by size and retrieve: length, total length, euclidean distance, total euclidean distance, branches, junctions, endpoints and number of filtered branches   

### Dependencies

* **poreTemplate.tif** file that should be downloaded from this repository and dropped into images directory
* **Multi-Template-Matching** + **IJ-OpenCV** Fiji plugins (installable via Update Sites)
* **3D ImageJ Suite**
* **IJPB** Fiji plugin (installable via Update Sites)

### Version history

Version 1 released on 08 january, 2026.

