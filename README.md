# Sentinel-1 Analysis Ready Data Preparation in Google Earth Engine
Google Earth Engine has been a valuable resource for earth observation researchers as it provided a big data platform for large area monitoring. Sentinel-1 image in-turn provided a global all weather day and night coverage which is suitable for large scale land monitoring. However, to preserve information content and user freedom, some value added processing steps are not applied on the images ingested in Google Earth Engine. In this technical note, we have implemented additional border noise correction, speckle filtering and radiometric terrain normalization to the Sentinel-1 data catalogue for land monitoring applications, that enables the data to be ready for immediate information extraction and interpretation.

![flowchart](https://user-images.githubusercontent.com/48068921/113488349-ce8cdd00-94bd-11eb-8bc5-a72f0fc7ecc6.png)

# Usage 
The Javascript implementation can be copied to GEE code editor. The python API provided here can be run using the terminal. The user should install earthengine-api locally or in the cloud VM. 

# Requirement
