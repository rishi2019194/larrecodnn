# CVN (Convolutional Visual Network) for PID in various LAr experiments

This is a common framework adapted from the original code developed for DUNE by Leigh W. and Saul A. M. See [*B. Abi et al. (DUNE Collaboration), "Neutrino interaction classification with a convolutional neural network in the DUNE far detector", Phys. Rev. D 102, 092003, 2020*](https://link.aps.org/doi/10.1103/PhysRevD.102.092003). The original codebase is situated [here](https://github.com/DUNE/dunereco/tree/develop/dunereco/CVN)

Key features of this framework include : 

1. A common producer class (`interfaces/PixelMapProducer` and `interfaces/ICVNMapper`) for creating wire-tick 2D images (pixel maps) for training based on a variety of inputs (hits, waveforms, true energy deposits). Example producer modules are defined in `modules/LArCVN*Mapper` 
2. Evaluation is done by the `modules/LArCVNEvaluator_module.cc` module that utilizes a tensorflow interfaces defined in `tools/TFNetHandler_tool.cc`
3. Dataset preparation for training is handled by the `modules/CVNZlibMaker_module.cc` module that processes the simulated input files and produces a compressed binary (`.gz`) of the pixel map per event as well as a `.info` file containing various truth information needed for training.
3. Various utility classes in `func/` are provided for preparing the training dataset as well as evaluation. Notably, `func/TrainingData.*` contains a template class for a user-defined output data structure that goes into the `.info` files. These include the training labels (flavor, interaction type etc) as well as custom auxillary information for general analysis purposes.  

Once the dataset is prepared for training, one can refer to (here)[https://github.com/DUNE/dunereco/tree/develop/dunereco/CVN/keras_scripts] for an example python framework that actually performs the training using the DUNE CVN architecture in keras-tensorflow for both GPU (highly recommended) and CPU-based setups.  
