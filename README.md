### Universality - SINDy
Test theories about universality in network dynamic with SINDy. The process is:
1. Load a graph.
  - Generate random graph (erdos-reyni or scale-free)
  - Load graph from published networks
2. Generate time-series data
  - epidemic ![epidemic](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bdx_i%7D%7Bdt%7D%3D-Bx_i%5Eb%2B%5Csum_%7Bj%3D1%7D%5E%7BN%7DA_i_jR%281-x_i%29x_j&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)
  - biochemical ![biochemical](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bdx_i%7D%7Bdt%7D%3DF-Bx_i-%5Csum_%7Bj%3D1%7D%5E%7BN%7DA_i_jRx_ix_j&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)
  - birth-dead ![birth-dead](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bdx_i%7D%7Bdt%7D%3D-Bx_i%5Eb%2B%5Csum_%7Bj%3D1%7D%5E%7BN%7DA_i_jRx_j%5Ea&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)
3. Fit SINDy model
