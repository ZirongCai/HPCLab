### Getting the source code
Clone this repository to your remote machine and switch to the project root directory:

$> cd 2dconv

### Environment
First of all, you need to set up an appropriate environment. The root directory contains a script called *conv_env.sh*.
Execute it as following:

$> source ./conv_env.sh

### Project Configuration for debugging:
Create a build directory:

$> mkdir build && cd build

Configure your project as following:

$> cmake ..

The default configuration allows to display images. Thus, it makes it is easy to debug the code. Ensure that you allowed
X11 forwarding to your ssh connection. Otherwise, you won't be able to see images on a screen of your local computer.
Type *make* in order to build your project:

$> make

### Project Configuration for performance testing:
You have to test your implementation on a production node once you are done with vectorization. Re-configure your
project as following:

$> cmake .. -DGRAPHICS=off

Build your project again:

$> make

Now, you have to submit a job to a production node. Fortunately, our cmake knows how to generate a proper job script for
you. If you've already configured the project then you will find a file **job.sh**. This script runs your code using
both 'default' and 'intrinsics' modes.  Execute it as following:

$> sbatch ./job.sh

Once a production node finishes executing of your task you will see a file called something **2Dconv.***. This file
contains your results. Open it and calculate your achieved speed-up.

