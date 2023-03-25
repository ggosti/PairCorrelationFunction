# PairCorrelationFunction
 Simple python code to compute the pair correlation function g(r) ([radial distribution function](https://en.wikipedia.org/wiki/Radial_distribution_function)).

Simple code that computes correctly and fast pair correlations.

## Usage

It requires an image with labels.

```python
centersX,centersY = find_centers(labels)
bs = 1000
insideBox,xsBox,ysBox=makeBox(centersX,centersY, bs)
rs = distances(centersX,centersY,xsBox,ysBox)
a = gdr(centersX,xsBox,rs,bs)
plt.plot(a[0], a[1])
plt.grid()
plt.hlines(1, 0,1000, color = "black", ls = "--")

```

# Other Projects

Here are some other projects for future comparrison.

https://github.com/aberut/paircorrelation2d

https://github.com/For-a-few-DPPs-more/structure-factor

https://github.com/mjo22/spatialstats

https://github.com/wenyan4work/point_cloud
