### Universality
1. Calculate local correlation matrix for biochemical network
2. Calculate correlation matrix using local matrix

### Universality Test
- Running multiple universality for normalization

### Universality - SINDy
Test theories about universality in network dynamic with SINDy. The process is:
1. Load a graph.
  - Generate random graph (erdos-reyni or scale-free)
  - Load graph from published networks
2. Generate time-series data
  - epidemic
  - biochemical
  - birth-dead
3. Fit SINDy model

### Information Flow -SINDy
Finding flow pattern in network dynamic with SINDy.
1. Generate time-series data for a known dynamic (biochemical)
2. Fit SINDy model with custom library
3. Find steady state in predicted network
  - Compare with real steady state
  - (**todo**) Compare estimated steady state followed by universality article
4. Apply a perturbation on steady state
5. Monitor information flow in network
6. (**todo**) Check universal pattern with calculated flow

### Seismology
Information flow in seismology data

### Financial
Predict behaivour of companies status in financial market

