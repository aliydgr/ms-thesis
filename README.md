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
