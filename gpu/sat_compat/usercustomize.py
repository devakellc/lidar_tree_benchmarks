# SegmentAnyTree (#M6) stack-modernization compat shims, applied as a
# `usercustomize` module: any python started with this directory on PYTHONPATH
# auto-imports it at interpreter startup (site machinery), so the upstream
# eval.py / torch_points3d run UNMODIFIED. The driver puts this dir on the
# PYTHONPATH it hands to run_inference.sh.
#
# WHY: SAT pins torch 1.9 / CUDA 11.1 / hydra 1.0 / omegaconf 2.0 / old PyG. The
# RTX 5090 (sm_120) forces torch 2.7 / cu128 / py3.11, which drags hydra 1.3 +
# omegaconf 2.3 + modern PyG. That version jump breaks several runtime
# assumptions. Each shim below restores ONE behavior SAT was written against.
#
# STATUS (2026-06-09): together with the image's pinned deps (torch-geometric
# 2.3.1, hdbscan, pyarrow) and --shm-size/--ipc=host, these shims carry SAT
# inference END-TO-END on the 5090 (verified: SOAP-native clip -> labeled .laz
# with PredInstance). The four shims below each restore one behavior the 2023-era
# stack assumed; the torch-geometric pin (not a shim) fixes the Data.keys drift.
import warnings


def _patch_omegaconf_missing_attr():
    # omegaconf 2.0.6: reading a MISSING attribute on a struct-disabled config
    # returns None (the "optional config group" idiom, e.g. self._cfg.training).
    # omegaconf >= 2.2 raises ConfigAttributeError instead. torch_points3d relies
    # on the old behavior pervasively; restore it in one place.
    try:
        import omegaconf
        from omegaconf import OmegaConf, DictConfig
    except Exception:
        return
    if getattr(DictConfig, "_sat_compat_getattr", False):
        return
    _orig = DictConfig.__getattr__

    def _getattr(self, key):
        try:
            return _orig(self, key)
        except omegaconf.errors.ConfigAttributeError:
            if not OmegaConf.is_struct(self):
                return None
            raise

    DictConfig.__getattr__ = _getattr
    DictConfig._sat_compat_getattr = True


def _patch_torch_load_weights_only():
    # torch 2.6 flipped torch.load's weights_only default to True, which refuses
    # the omegaconf run_config pickled into the SAT checkpoint. The official
    # weights are trusted -> restore weights_only=False.
    try:
        import torch
    except Exception:
        return
    if getattr(torch, "_sat_compat_load", False):
        return
    _orig = torch.load

    def _load(*a, **k):
        k.setdefault("weights_only", False)
        return _orig(*a, **k)

    torch.load = _load
    torch._sat_compat_load = True


def _set_spawn_start_method():
    # SAT's meanshift/hdbscan clustering builds a multiprocessing.Pool AFTER the
    # model forward has initialized CUDA in the parent. The Linux default 'fork'
    # inherits the live CUDA context + torch's thread pool into the child, which
    # deadlocks (idle parent, zombie workers). 'spawn' starts a clean
    # interpreter -> CUDA-safe. SAT runs num_workers=0, so no DataLoader workers
    # are affected. Set at interpreter startup, before any Pool is created.
    try:
        import multiprocessing as mp
        if mp.get_start_method(allow_none=True) != "spawn":
            mp.set_start_method("spawn", force=True)
    except Exception:
        pass


def _patch_tensor_getitem_cross_device():
    # torch < ~2.1 auto-moved a device-mismatched index tensor; torch 2.7 raises
    # "indices should be ... on the same device as the indexed tensor". SAT's
    # 2023 code indexes CPU tensors with GPU index/mask tensors (and vice versa)
    # in several places (region_grow, PointGroup _compute_score, ...). Restore
    # the old leniency once: when a Tensor key's device != the tensor's, move the
    # key to the tensor's device, including tuple/list indexing forms such as
    # x[idx, :] and x[:, mask]. Fast path (non-Tensor keys) is untouched.
    try:
        import torch
    except Exception:
        return
    if getattr(torch.Tensor, "_sat_compat_getitem", False):
        return
    _orig = torch.Tensor.__getitem__

    def _normalize_key(key, device):
        if isinstance(key, torch.Tensor) and key.device != device:
            return key.to(device)
        if isinstance(key, tuple):
            return tuple(_normalize_key(k, device) for k in key)
        if isinstance(key, list):
            return [_normalize_key(k, device) for k in key]
        return key

    def _getitem(self, key):
        return _orig(self, _normalize_key(key, self.device))

    torch.Tensor.__getitem__ = _getitem
    torch.Tensor._sat_compat_getitem = True


def apply():
    _patch_omegaconf_missing_attr()
    _patch_torch_load_weights_only()
    _set_spawn_start_method()
    _patch_tensor_getitem_cross_device()


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    apply()
