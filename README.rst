=======================
Localrec scipion plugin
=======================

Electron cryomicroscopy can yield near-atomic resolution structures of highly ordered macromolecular complexes. Often however some subunits bind in a flexible manner, have different symmetry from the rest of the complex, or are present in sub-stoichiometric amounts, limiting the attainable resolution. Here we report a general method for the localized three-dimensional reconstruction of such subunits. (see http://www.opic.ox.ac.uk/wiki/index.php/Localized_reconstruction for details) 

This sub-package contains data and protocol classes to use Localrec within the Scipion framework


=====
Setup
=====

- **Install this plugin:**

.. code-block::

    scipion installp -p scipion-em-localrec

OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

Alternatively, in devel mode:

.. code-block::

    scipion installp -p local/path/to/scipion-em-localrec --devel
