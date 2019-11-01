_CORES?=$(shell grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
_PREFIX?=?/home/blaze/install
_MALLOC?=system
_BUILD_TYPE?=Debug
_DOCKER_FILE?=Dockerfile
_SRC?=/home/balze/bblaze/

build:
	docker build -t bblaze                \
	--build-arg CORES=$(_CORES)  	      \
	--build-arg MALLOC=$(_MALLOC)	      \
	--build-arg PREFIX=$(_PREFIX)	      \
	--build-arg BUILD_TYPE=$(_BUILD_TYPE) \
	-f $(_DOCKER_FILE)

