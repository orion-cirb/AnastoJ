# ProtruJ_Morphometrics

* **Developed by:** Thomas & Héloïse
* **Developed for:** Laurent & Emmanuelle
* **Team:** Germain
* **Date:** May 2025
* **Software:** Fiji


### Images description

3D images taken with x10 objective on a spinning-disk microscope.

3 channels:
  1. DAPI nuclei (not used)
  2. GFP protrusions
  3. Brightfield pores

### Macro description

* Ask the user to rotate the image so the pores are pointing to the right
* Detect pores in the 3rd channel using Multi-Template-Matching
* Segment protrusions in 2nd channel
* For each detected pore:
     * Crop the protrusions binary mask around it, filter out small objects, and clear the interior of the pore
     * Run skeletonization + geodesic distance map + local thickness from the resulting mask
     * For each protrusion (= skeleton), compute the following: nearest corner of the pore, maximum distance from the pore, direction of the centroid, and mean thickness

 #### Pore corners positioning

<img src="https://github.com/user-attachments/assets/9ee14349-463f-43ec-bea9-94cda6894ec1" width="200">


### Dependencies

* **poreTemplate.tif** file that should be downloaded from this repository and dropped into images directory
* **Multi-Template-Matching** + **IJ-OpenCV** Fiji plugins (installable via Update Sites)
* **IJPB** Fiji plugin (installable via Update Sites)

### Version history

Version 1 released on May 14, 2025.

