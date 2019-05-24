=======================
Localrec scipion plugin
=======================

Electron cryomicroscopy can yield near-atomic resolution structures of highly ordered macromolecular complexes. Often however some subunits bind in a flexible manner, have different symmetry from the rest of the complex, or are present in sub-stoichiometric amounts, limiting the attainable resolution. Here we implement a general method for the localized three-dimensional reconstruction of such subunits. (see http://www.opic.ox.ac.uk/wiki/index.php/Localized_reconstruction for details) 

This sub-package contains data and protocol classes to use Localrec within the Scipion framework


===================
Install this plugin
===================

You will need to use `2.0.0 <https://github.com/I2PC/scipion/releases/tag/v2.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block::

    scipion installp -p scipion-em-localrec

OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-localrec.git

2. Install:

.. code-block::

           scipion installp -p path_to_scipion-em-localrec --devel

- **Tests**

1. scipion test localrec.tests.test_protocol_localized_reconstruction

========
Protocols
========

* localized subparticles: Calculate the orientations of the subunits of interest and their positions in the original particle images.
* filter_subunits: Filter the subunits (sub-particles) based on spatial distance, angular distance, etc.
* localized extraction: Extract computed sub-particles from a SetOfParticles.


===============
Buildbot status
===============
TODO: add to buildbot
Status devel version: 

.. image:: http://arquimedes.cnb.csic.es:9980/badges/ccp4_devel.svg

Status production version: 

.. image:: http://arquimedes.cnb.csic.es:9980/badges/ccp4_prod.svg




