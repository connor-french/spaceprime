# Suitability to deme size transformation functions
## Brief overview
All transformations are performed on a raster with 0-1 suitability values (or a raster that will be normalized to 0-1). The output of the transformation is an array of the same shape as the input raster, but with values transformed according to the specified function.

The "linear" transformation multiplies the input values by the maximum local deme size (`max_local_size`).

The "threshold" transformation creates an array where values below the `threshold` value are set to zero and values about the `threshold` value are set to 1.

The "sigmoid" transformation applies a sigmoid function to the data using Eq. 1 from [Frazier and Wang 2013, Modeling landscape structure response across a gradient of land cover intensity](https://www.researchgate.net/publication/257319938_Modeling_landscape_structure_response_across_a_gradient_of_land_cover_intensity), where an `inflection_point` and `slope` are specified. The `inflection_point` can be thought of like a `threshold` value, where original values below this value descend quicker to zero, and values about this value increase quicker to 1. The slope determines how fast values change on either side of the inflection point. A sufficiently steep slope makes this a threshold function, while a sufficiently shallow slope makes this a linear function. 

