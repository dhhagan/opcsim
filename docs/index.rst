.. raw:: html

    <style type="text/css">
    .thumbnail {{
        position: relative;
        float: left;
        margin: 10px;
        width: 180px;
        height: 200px;
    }}

    .thumbnail img {{
        position: absolute;
        display: inline;
        left: 0;
        width: 170px;
        height: 170px;
    }}

    </style>

opcsim: Simulating Optical Particle Sensors
============================================

.. raw:: html

    <div style="clear: both"></div>

    <div class="container-fluid hidden-xs hidden-sm">
      <div class="row">
        <a href="examples/three_weights.html">
          <div class="col-md-2 thumbnail">
            <img src="_static/three_weights_thumb.png">
          </div>
        </a>
        <a href="examples/urban_distribution_pdf.html">
          <div class="col-md-2 thumbnail">
            <img src="_static/urban_distribution_pdf_thumb.png">
          </div>
        </a>
        <a href="examples/ten_bin_opc.html">
          <div class="col-md-2 thumbnail">
            <img src="_static/ten_bin_opc_thumb.png">
          </div>
        </a>
        <a href="examples/opc_with_dist_number_and_vol.html">
          <div class="col-md-2 thumbnail">
            <img src="_static/opc_with_dist_number_and_vol_thumb.png">
          </div>
        </a>
        <a href="examples/opc_with_dist.html">
          <div class="col-md-2 thumbnail">
            <img src="_static/opc_with_dist_thumb.png">
          </div>
        </a>
        <a href="examples/build_your_own_distribution.html">
          <div class="col-md-2 thumbnail">
            <img src="_static/build_your_own_distribution_thumb.png">
          </div>
        </a>
      </div>
    </div><!-- /.container-fluid -->

    <br>

    <div class="container-fluid">
      <div class="row">
        <div class="col-md-6">
        
          <br>

`opcsim` is a Python library for simulating the response of low-cost optical
particle sensors (OPS's) to various aerosol distributions to better understand
the sources of error and limitations of these devices. It provides an
easy-to-use API for building simple OPC and Nephelometer models as well as 
model and visualize aerosol distributions.

For more information, please read our paper available in `Atmospheric Measurement Techniques <https://coming-soon-to-amt.com>`_.


To view the source code or report a bug, please visit the `github repository
<https://github.com/dhhagan/opcsim>`_.


.. raw:: html

   </div><!-- /.col-md-6 -->

   <div class="col-md-3">
    <div class="panel panel-default">
      <div class="panel-heading">
        <h3 class="panel-title">Contents</h3>
      </div>

      <div class="panel-body">

.. toctree::
   :maxdepth: 1

   Installation <installing>
   Contributing <contributing>
   API reference <api>
   Tutorial <tutorial>
   Example gallery <examples/index>

.. raw:: html

      </div><!-- /.panel-body -->
    </div><!-- /.panel panel-default -->
  </div><!-- /.col-md-3 -->

  <div class="col-md-3">
    <div class="panel panel-default">
      <div class="panel-heading">
        <h3 class="panel-title">Features</h3>
      </div>
  
      <div class="panel-body">

* Simulate optical particle counters
* Simulate nephelometers
* Understand how particle sensors react to changes in aerosol size and composition
* Easily visualize aerosol distirbutions

.. raw:: html

      </div><!-- /.panel-body-->
    </div><!-- /.panel-->
   </div><!-- /.col-md-3-->

  </div><!-- /.row -->

  <div class="row">


.. raw:: html

    <div class="col-md-9">

    <h3>Abstract</h3>

    <p>
    Low-cost sensors for measuring particulate matter (PM) offer the ability to understand
    human exposure to air pollution at spatiotemporal scales that have previously been
    impractical. However, such low-cost PM sensors tend to be poorly characterized, and
    their measurements of mass concentration can be subject to considerable error. Recent
    studies have investigated how individual factors can contribute to this error, but these
    studies are largely based on empirical comparisons and generally do not examine the
    role of multiple factors simultaneously. Here, we present a new physics-based framework
    and open-source software package (opcsim) for evaluating the ability of low-cost optical
    particle sensors (optical particle counters and nephelometers) to accurately characterize
    the size distribution and/or mass loading of aerosol particles. This framework, which uses
    Mie Theory to calculate the response of a given sensor to a given particle population, is
    used to estimate the relative error in mass loading for different sensor types, given
    variations in relative humidity, aerosol optical properties, and the underlying particle size
    distribution. Results indicate that such error, which can be substantial, is dependent on
    the sensor technology (nephelometer vs. optical particle counter), the specific
    parameters of the individual sensor, and differences between the aerosol used to
    calibrate the sensor and the aerosol being measured. We conclude with a summary of
    likely sources of error for different sensor types, environmental conditions, and particle
    classes, and offer general recommendations for choice of calibrant under different
    measurement scenarios.
    </p>

    <br>

    To cite this work, please use the following:

    <h4 style="line-height: 1.2em; margin-bottom:25px;">
    Hagan, D.H. and Kroll, Jesse H.: Assessing the accuracy of low-cost optical particle sensors using a 
    physics-based approach, Atmos. Meas. Tech. Disc., submitted, 2020.
    </h4>

    </div><!-- /.col-md-9-->
  
  </div><!-- /.row -->
  </div><!-- /.container-fluid -->

