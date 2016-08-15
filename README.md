# Piranhas
A Toolkit for creating Peripheral Architectures.

![FoveatedRepresentations](/images/Total_Foveated_Cartoon.png)

This toolkit was implemented and used to create the peripheral architectures of the papers:
* **_Can Peripheral Representations Improve Clutter Metrics on Complex Scenes?_**. Arturo Deza & Miguel P. Eckstein, NIPS 2016.
* **_Object Detection Through Exploration With A Foveated Visual Field_**. Emre Akbas & Miguel P. Eckstein, ArXiv 2014 (under review).

The toolkit was originally written in MATLAB, but has been extended to python and Torch to increase cross-collaborations between fields of vision science, computer vision and deep learning. Piranhas stands for Peripheral Architectures for Natural, Hybrid and Artificial Systems, we decided to make our toolkit public to stimulate possible 'hybrid' ideas in the general vision (human,computer,robot) community.

# What is a Peripheral Architecture?
A peripheral architecture is a collection of regions that simulate human-like pooling regions and foveal and peripheral like mechanisms. Peripheral Architectures can be simple such as the ones used in "Multiple Object Recognition with Visual Attention" [Ba,Mnih & Kavukucoglu, 2015], or can be more complex given biological constraints such as the one proposed in "Metamers of the Ventral Stream" [Freeman & Simoncelli, 2011].

![Observer1](http://imgur.com/uhw8Lq3.gif) ![Observer2](http://imgur.com/SL0WTgH.gif) ![Observer3](http://imgur.com/SL0WTgH.gif)
![Observer4](http://imgur.com/uhw8Lq3.gif) ![Observer5](http://imgur.com/SL0WTgH.gif) ![Observer6](http://imgur.com/SL0WTgH.gif)

A sample collection of 6 human observers doing a single trial of a target search task. Real EyeTracking Data on the right of each subfigure, and on the left Foveated Feature Congestion(FFC) (Deza & Eckstein, 2016) is superimposed given every fixation. Similar visualizations can be extended to computer generated fixations or computer vision models with foveated systems in search (Akbas & Eckstein, 2014). Visualizations were generated using the Piranhas toolkit.

# Using Piranhas

1. Download the Piranhas toolbox to create Peripheral Architectures at your convenience in [Matlab](https://github.com/ArturoDeza/Piranhas/tree/master/MATLAB), [python](https://github.com/ArturoDeza/Piranhas/tree/master/python) or [Torch](https://github.com/ArturoDeza/Piranhas/tree/master/torch):

2. Define your Computer + Human perception parameters. Read the [Tutorial](https://github.com/ArturoDeza/Piranhas/tree/master/Tutorial) to learn about these parameters.

3. Create a Peripheral Architecture.

4. Pool your dense feature maps (Currently availabe in MATLAB, Coming soon for python and torch)

## [FAQ]: Frequenty Asked Questions:

### Q: What are examples of papers that mix (computer/human/robot) perception?

A: Some papers that have a mixed human + computer perception flavor, be it through experiments or discussion are:

* Agrawal, Pulkit, et al. "Pixels to voxels: modeling visual representation in the human brain." arXiv preprint arXiv:1407.5104 (2014).
* Aminoff, Elissa M., et al. "Applying artificial vision models to human scene understanding." Frontiers in computational neuroscience 9 (2015).
* Pramod, R. T., and S. P. Arun. "Do computational models differ systematically from human object perception?." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2016.
* Das, Abhishek, et al. "Human Attention in Visual Question Answering: Do Humans and Deep Networks Look at the Same Regions?." arXiv preprint arXiv:1606.03556 (2016).
* Kriegeskorte, Nikolaus. "Deep neural networks: A new framework for modeling biological vision and brain information processing." Annual Review of Vision Science 1 (2015): 417-446.
* Cichy RM, Khosla A, Pantazis D, Torralba A, Oliva A (2016). Comparison of deep neural networks to spatio-temporal cortical dynamics of human visual object recognition reveals hierarchical correspondence SciReports, 6:27755. doi: 10.1038/srep27755.
* Can Peripheral Representations Improve Clutter Metrics on Complex Scenes? Arturo Deza & Miguel P. Eckstein, NIPS 2016.
* Object Detection Through Exploration With A Foveated Visual Field. Emre Akbas & Miguel P. Eckstein, (under review).

Do you have a paper at the intersection of more than one vision field? Let us know to expand this preliminary list!

### Q: I want to create my own peripheral architecture, what should I do?

Once you have a functional architecture, you are more than welcome to add it to the Piranha School (similar to the Caffe Model Zoo). 
Coming soon!

### Q: I do Deep Learning, what can I (get out of) / (do for)  Piranhas toolkit?

Deep Learning has begun embracing a strong trend in mimicking neuroscience-like mechanisms in their algorithms. Take for example: Perceptrons with Neurons, hierarchical structure of V1, V2, V4, and IT with Deep Networks with CNN's, and RNN's + LSTM's with Overt Attention Mechanisms.
While there is still much more work to be done in Deep Learning, and we are far from understanding intelligence, there have been remarkable successes
in borrowing these ideas from Neuroscience and giving them an engineering twist to make them applicable in Deep Learning systems.

We believe that one of the main unexplored areas in Computer Vision and Deep Learning is the use of a periphery and/or covert systems in recognition, registration and possibly egocentric vision. Perhaps with this toolkit, many deep learning scientists and engineers can help the vision community understand how/why/if peripheral architectures are useful in computer vision.

### Q: I do Computer Vision, what can I (get out of) / (do for) Piranhas toolkit?

One of the main goals in Computer Vision is to give the gift of sight to computers i.e. produce human-like levels of performance in ILSVRC and on image datasets like MNIST, LabelMe, SUN, ImageNet, and MSCOCO. While it is debatable if computers should in fact have the same (artificial) visual perception algorithms as humans, we think that the use of experimenting with cleverly designed peripheral architectures opens a realm of possibilities for future solutions and new problems in the field of computer vision. 

### Q: I do Robotics, what can I (get out of) / (do for) Piranhas toolkit?

Robotics is a curious field, where Machine Vision plays a critical role in the general pipeline that includes developing algorithms for locomotion, planning, and sensing. Creating robots with human or animal like characteristics is a current hot topic, and adding peripheral constraints similar to an animals' eye may play an important role in path-planning for an agent in its environment, and may lead to intuitions of how an animal might explore an environment given its (artificially created) visual limitations. Fun fact: Did you that the hummingbird has two [foveas](https://en.wikipedia.org/wiki/Fovea_centralis)?

# Credits & Citation
Piranhas Toolkit was mainly written by Arturo Deza, and Emre Akbas. All code was written inside the 
[VIU-Lab @ UCSB](https://labs.psych.ucsb.edu/eckstein/miguel/index.htm).

If you found this code useful for your research, please cite:
```
@misc{Deza2016piranhas,
author = {Deza, Arturo and Abkas, Emre and Eckstein, Miguel P.},
title = {Piranhas Toolkit: Peripheral Architectures for Natural, Hybrid and Artificial Systems}
year = {2016},
publisher = {GitHub},
journal = {GitHub repository},
howpublished = {\url{https://github.com/ArturoDeza/Piranhas}}
}
```

Piranhas is released under the [BSD3-license](https://github.com/ArturoDeza/Piranhas/blob/master/license.txt).

README.md file written by Arturo Deza.
