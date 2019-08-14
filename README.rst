=======================
Localrec Scipion plugin
=======================

Cryogenic electron microscopy (cryoEM) can yield near-atomic resolution structures of highly ordered macromolecular complexes. Often however some subunits bind in a flexible manner, have different symmetry from the rest of the complex, or are present in sub-stoichiometric amounts, limiting the attainable resolution. We have developed a general method for the localized three-dimensional reconstruction of such subunits earlier (Ilca et al 2015 Nature Commun).

In this project we implement Localized Reconstruction as a plugin (Locarec) for Scipion.

If Localized Reconstruction is useful in your work, please cite:

Ilca SL, Kotecha A, Sun X, Poranen MM, Stuart DI & Huiskonen JT (2015). Localized reconstruction of subunits from electron cryomicroscopy images of macromolecular complexes. Nat Commun 6, 8843. doi:10.1038/ncomms9843


===================
Install this plugin
===================

You will need to use `2.0.0 <https://github.com/I2PC/scipion/releases/tag/v2.0>`_ version of Scipion or later to run these protocols. Two versions of the plugin are available:

**Stable version**  

Either use the following command

.. code-block::

    scipion installp -p scipion-em-localrec

or

Install the plugin through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

**Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-localrec.git

2. Install:

.. code-block::

           scipion installp -p path_to_scipion-em-localrec --devel

3. Test the plugin:

.. code-block::

           scipion test localrec.tests.test_protocol_localized_reconstruction

========
Protocols
========

The following protocols have been implemented so far

* localrec - define subparticles: Calculate the orientations of the sub-particlesunits of interest and their positions in the original particle images.
* localrec - filter subparticles: Filter the sub-particles based on spatial distance, angular distance, etc.
* localrec - extract subparticles: Extract computed sub-particles from a SetOfParticles.
* localrec - stitch subparticles: Regenerate a particle (volume) from its sub-particle based reconstructions. 

===============
Buildbot status
===============
Status devel version: 

.. image:: http://arquimedes.cnb.csic.es:9980/badges/ccp4_devel.svg

Status production version: 

.. image:: http://arquimedes.cnb.csic.es:9980/badges/ccp4_prod.svg


================
Acknowledgements
================

The authors acknowledge funding from the European Research Council under the European Union’s Horizon 2020 research and innovation programme (649053).
