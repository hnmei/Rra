<!--
 * @Author: Haonan Mei
 * @Date: 2022-06-03 15:29:20
 * @LastEditTime: 2022-06-06 23:28:12
 * @LastEditors: Haonan Mei
 * @Description: 
 * @FilePath: \undefinedd:\leslie\ra\taozeng\00_rra_to_r\working\Rra\README.md
-->

# Rra
R version implementation of the reduced-rank approach by Dashan Huang et al. (MS, 2022)

## Install
```{R}
library(devtools)
devtools::install_github("hnmei/Rra")
```

## Data
We provide two data sets:

**"48_ports_ret_74-16.Rdata":**
48 industry portfolios (516\*49 dataframe, the added one is date) as the example input of return data.

**"70_factor_proxies_74-16.Rdata":**
70 factor proxies (516\*71 dataframe, the added one is date) as the example input of factors.

**`date` column is needed and if you use your own data, notice that put the `date` column at the first column.**

## Usage
```{R}
# load dependent pkg
library(expm)
library(POET)
```

```{R}
library(Rra)

data("48_ports_ret_74-16") # return data of basis portfolios
data("70_factor_proxies_74-16") # factor proxies
rra <- rra(ff48, factor_proxies, 5)
```

### RRA with price restriction
The dafault price error(alpha) per year is 1%, adjusted with standard error of each portfolio automatically.
```{R}
rra_pere <- rra(ff48, factor_proxies, 5, price_error=TRUE)
```

Also, you can define alpha yourself, like with 5% price error per year:
```{R}
rra_pere <- rra(
    ff48, factor_proxies, 5, 
    price_error=TRUE, alpha=5
)
```

### RRA with changing weight
We estimate the weight matrix (S2) with the method of with the principal orthogonal complement thresholding method of Fan, Liao, and Mincheva (2013). The R code implementation can refer to R package `POET`. The weight matrix is estimated automatically.
```{R}
rra_w <- rra(
    ff48, factor_proxies, 5
    weight=TRUE
)
```

## References
1. Dashan Huang et al. (MS 2022) Shrinking Factor Dimension A Reduced-Rank Approach
2. Fan, J., Liao, Y., Mincheva, M., 2013. Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society. Series B (Statistical Methodology) 75, 603–680.


