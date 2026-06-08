ARG PYTORCH=2.7.0
ARG CUDA=12.8
ARG CUDNN=9

FROM pytorch/pytorch:${PYTORCH}-cuda${CUDA}-cudnn${CUDNN}-devel

ENV TORCH_NVCC_FLAGS="-Xfatbin -compress-all"
ARG CUDA_SM_LIST="6.0 6.1 6.2 7.0 7.2 7.5 8.0 8.6 8.9 9.0 12.0"
ENV CMAKE_CUDA_ARCHITECTURES=${CUDA_SM_LIST} \
    TORCH_CUDA_ARCH_LIST=${CUDA_SM_LIST}

RUN apt-get update && apt-get install -y --no-install-recommends \
    git ninja-build cmake build-essential libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
RUN git clone --depth 1 --branch cuda12-installation \
    https://github.com/CiSong10/MinkowskiEngine.git

WORKDIR /workspace/MinkowskiEngine

# Fix std::to_address ambiguity (issue #596); path varies by base image gcc.
RUN for f in /usr/include/c++/11/bits/shared_ptr_base.h /usr/include/c++/12/bits/shared_ptr_base.h; do \
      if [ -f "$f" ]; then \
        sed -i 's/\bauto __raw = __to_address(__r.get());/auto __raw = std::__to_address(__r.get());/' "$f"; \
      fi; \
    done

RUN python setup.py install --force_cuda --blas=openblas

COPY test_me.py /workspace/test_me.py
CMD ["python", "/workspace/test_me.py"]
