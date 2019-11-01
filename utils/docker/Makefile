_CORES?=$(shell nproc)
_PREFIX?=/home/blaze/install
_MALLOC?=system
_BUILD_TYPE?=Debug
_DOCKER_FILE?=Dockerfile
_SRC?=/home/blaze/bblaze/

build:
	docker build -t batched_blaze         \
	--build-arg SRC=$(_SRC)  	          \
	--build-arg CORES=$(_CORES)  	      \
	--build-arg MALLOC=$(_MALLOC)	      \
	--build-arg PREFIX=$(_PREFIX)	      \
	--build-arg BUILD_TYPE=$(_BUILD_TYPE) \
	-f $(_DOCKER_FILE)                    \
	.

force-build:
	docker build --no-cache -t batched_blaze \
	--build-arg SRC=$(_SRC)  	          	 \
	--build-arg CORES=$(_CORES)  	      	 \
	--build-arg MALLOC=$(_MALLOC)	      	 \
	--build-arg PREFIX=$(_PREFIX)	      	 \
	--build-arg BUILD_TYPE=$(_BUILD_TYPE) 	 \
	-f $(_DOCKER_FILE)                    	 \
	.
