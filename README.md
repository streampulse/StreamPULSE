# **StreamPULSE**

### An R package for fitting metabolism models to StreamPULSE data
### **Description**
StreamPULSE is a community of researchers working to
answer fundamental questions about ecosystem production and
respiration in streams and rivers. The StreamPULSE project
[streampulse.org](http://www.streampulse.org)
includes online tools for viewing, archiving, and cleaning data
submitted by researchers from around the globe. It leverages 
metabolism modeling frameworks
[streamMetabolizer](https://github.com/USGS-R/streamMetabolizer)
and [BASE](https://github.com/dgiling/BASE).
This package facilitates 
data acquisition from the StreamPULSE database, data
supplementation where necessary, formatting for streamMetabolizer
or BASE, model fitting and prediction, and visualization of model outputs.

---
### **Installation**
```R
# in R
library(devtools)
install_github('streampulse/StreamPULSE')
```

---
### **Example**
See StreamPULSE's [https://data.streampulse.org/model](model page)
for an example of how to use this package.
