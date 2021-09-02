# Installation (UBUNTU)

## Libraries

### GLFW (Graphics Library Framework)

https://stackoverflow.com/questions/17768008/how-to-build-install-glfw-3-and-use-it-in-a-linux-project

Download the latest version (source package) from http://www.glfw.org/download.html and unzip anywhere on your system followed by

```console
foo@bar:~$ cmake -G "Unix Makefiles"
foo@bar:~$ make
foo@bar:~$ sudo make install
```

This will install the library into `usr/local/include/GLFW` and `/usr/local/lib/` which are most likely already in your paths file.

### GLM (OpenGL Mathematics)

Download the latest version from https://github.com/g-truc/glm/tags and unzip anywhere on your system followed by

```console
foo@bar:~$ cmake -G "Unix Makefiles"
foo@bar:~$ make
foo@bar:~$ sudo make install
```

This will install the library into `/usr/local/include/glm` and `/usr/local/lib/` which are most likely already in your paths file.

If you get the following compilation error message
```console
/usr/local/include/glm/gtx/transform.hpp:23:3: error: #error "GLM: GLM_GTX_transform is an experimental extension and may change in the future. Use #define GLM_ENABLE_EXPERIMENTAL before including it, if you really want to use it."
```
then you need to define `GLM_ENABLE_EXPERIMENTAL` in the file `/usr/local/include/glm/detail/setup.hpp` by adding the line
```
#define GLM_ENABLE_EXPERIMENTAL 0
```
if you would like to turn off any expermental fetures and otherwise
```
#define GLM_ENABLE_EXPERIMENTAL 1
```

### GLAD (GL/GLES/EGL/GLX/WGL Loader-Generator)

Visit https://glad.dav1d.de/ to generate a loader.

To check the OpenGL (gl) version:
```console
foo@bar:~$ glxinfo | grep "OpenGL version"
```

Add useful extensions and leave the rest of the form as is.
Download the generated zip file and unzip in `libs/glad`.


### FreeType (FreeType is a freely available software library to render fonts.)

Download the latest version from https://freetype.org and unzip anywhere on your system followed by

```console
foo@bar:~$ cmake -E make_directory build
foo@bar:~$ cmake -E chdir build cmake ..
foo@bar:~$ make
foo@bar:~$ sudo make install
```

This will install the library into `/usr/local/include/freetype2` and `/usr/local/lib/` which are most likely already in your paths file.z



## Compile and setup libraries

Simply run `scripts/build_libs.sh`





















# Installation (MAC)

## Libraries

### GLFW (Graphics Library Framework)

https://stackoverflow.com/questions/17768008/how-to-build-install-glfw-3-and-use-it-in-a-linux-project

Download the latest version (source package) from http://www.glfw.org/download.html and unzip anywhere on your system followed by

```console
foo@bar:~$ cmake -G "Unix Makefiles"
foo@bar:~$ make
foo@bar:~$ sudo make install
```

This will install the library into `usr/local/include/GLFW` and `/usr/local/lib/` which are most likely already in your paths file.

### GLM (OpenGL Mathematics)

Download the latest version from https://github.com/g-truc/glm/tags and unzip anywhere on your system followed by

```console
foo@bar:~$ cmake -G "Unix Makefiles"
foo@bar:~$ make
foo@bar:~$ sudo make install
```

This will install the library into `/usr/local/include/glm` and `/usr/local/lib/` which are most likely already in your paths file.

If you get the following compilation error message
```console
/usr/local/include/glm/gtx/transform.hpp:23:3: error: #error "GLM: GLM_GTX_transform is an experimental extension and may change in the future. Use #define GLM_ENABLE_EXPERIMENTAL before including it, if you really want to use it."
```
then you need to define `GLM_ENABLE_EXPERIMENTAL` in the file `/usr/local/include/glm/detail/setup.hpp` by adding the line
```
#define GLM_ENABLE_EXPERIMENTAL 0
```
if you would like to turn off any expermental fetures and otherwise
```
#define GLM_ENABLE_EXPERIMENTAL 1
```

### GLAD (GL/GLES/EGL/GLX/WGL Loader-Generator)

Visit https://glad.dav1d.de/ to generate a loader.

To check the OpenGL (gl) version visit https://support.apple.com/en-us/HT202823

Add useful extensions and leave the rest of the form as is.
Download the generated zip file and unzip in `libs/glad`.


### FreeType (FreeType is a freely available software library to render fonts.)

Download the latest version from https://freetype.org and unzip anywhere on your system followed by

```console
foo@bar:~$ cmake -E make_directory build
foo@bar:~$ cmake -E chdir build cmake ..
foo@bar:~$ make
foo@bar:~$ sudo make install
```

This will install the library into `/usr/local/include/freetype2` and `/usr/local/lib/` which are most likely already in your paths file.



## Compile and setup libraries

Simply run `scripts/build_libs.sh`