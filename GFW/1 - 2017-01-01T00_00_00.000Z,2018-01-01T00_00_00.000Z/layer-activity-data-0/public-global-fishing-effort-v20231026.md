
## AIS fishing bq dataset
### Description
Global Fishing Watch uses data about a vessel’s identity, type, location, speed, direction and more that is broadcast using the Automatic Identification System (AIS) and collected via satellites and terrestrial receivers. AIS was developed for safety/collision-avoidance. Global Fishing Watch analyzes AIS data collected from vessels that our research has identified as known or possible commercial fishing vessels, and applies a fishing presence algorithm to determine “apparent fishing activity” based on changes in vessel speed and direction. The algorithm classifies each AIS broadcast data point for these vessels as either apparently fishing or not fishing and shows the former on the Global Fishing Watch fishing activity heat map. AIS data as broadcast may vary in completeness, accuracy and quality. Also, data collection by satellite or terrestrial receivers may introduce errors through missing or inaccurate data. Global Fishing Watch’s fishing presence algorithm is a best effort mathematically to identify “apparent fishing activity.” As a result, it is possible that some fishing activity is not identified as such by Global Fishing Watch; conversely, Global Fishing Watch may show apparent fishing activity where fishing is not actually taking place. For these reasons, Global Fishing Watch qualifies designations of vessel fishing activity, including synonyms of the term “fishing activity,” such as “fishing” or “fishing effort,” as “apparent,” rather than certain. Any/all Global Fishing Watch information about “apparent fishing activity” should be considered an estimate and must be relied upon solely at your own risk. Global Fishing Watch is taking steps to make sure fishing activity designations are as accurate as possible. Global Fishing Watch fishing presence algorithms are developed and tested using actual fishing event data collected by observers, combined with expert analysis of vessel movement data resulting in the manual classification of thousands of known fishing events. Global Fishing Watch also collaborates extensively with academic researchers through our research program to share fishing activity classification data and automated classification techniques.

Filters: geartype IN ('trawlers') and timestamp >= '2017-01-01T00:00:00.000Z' and timestamp <= '2018-01-01T00:00:00.000Z'
Resolution: low (10th degree resolution)

### Format
- GeoTIFF
- CRS: EPSG:4326 - WGS 84 - Geographic

### License
Unless otherwise stated, Global Fishing Watch data is licensed under a [Creative Commons Attribution-ShareALike 4.0 International license](https://creativecommons.org/licenses/by-sa/4.0/) and code under an [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0).
	