# PyColormap4Matlab
Simple Matlab and python script that import colormaps from matplotlib into Matlab.

*getPyPlot_cMap* returns any colormap implemented in the matplotlib python library. It calls a python script that writes the colormap matrix into a temporary file, i.e. **python** (and the matplotlib module) is **required**.
However, the advantage is that you get all the colormaps implemented in *matplotlib* and that you can specify the number of RGB quantization levels, i.e. the number of colors of the colormap.

A list of colormap names is provided in the function help section. `getPyPlot_cMap('!GetNames')` returns a cellstring containing all available colormap names.
See https://matplotlib.org/examples/color/colormaps_reference.html for an illustration of colormaps.

Also available on MatlabCentral File Exchange:  
https://de.mathworks.com/matlabcentral/fileexchange/68239-pycolormap4matlab


## Usage
`cMapNames = getPyPlot_cMap('!GetNames')`  
Returns a cellstring containing all available colormap names.

`cMap = getPyPlot_cMap(cMapName)`  
Returns the colormap *cMapName* with the default of 128 colors. cMap will be a 128x3 matrix.

`cMap = getPyPlot_cMap(cMapName, NumberOfColors)`  
Specify the number of colors, i.e. the number of rows in cMap.

`cMap = getPyPlot_cMap(cMapName, NumberOfColors, keepAlphaChannel)`  
If keepAlphaChannel is not 0 cMap has a 4th column containing the alpha channel.

`cMap = getPyPlot_cMap(cMapName, NumberOfColors, keepAlphaChannel, pythonSystemCommand)`  
Lets you specify the python command (possibly including a path, see below) used to execute the python script.


## Errors
If you have python installed but Matlab says `There was an error executing the command... System returned:...` you can try to pass the path to your python installation explicitly as the 4th parameter, e.g.:
```
cmp = getPyPlot_cMap('Accent', [], [], '"c:\Program Files\Python37\python.exe"');
```
Note the double-quotes around the path, which are neccessary because of the containing space character.
