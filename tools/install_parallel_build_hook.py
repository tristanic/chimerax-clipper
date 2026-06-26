"""Install (or remove) the parallel-build startup hook in ChimeraX's user site.

Run by ChimeraX so it can resolve its own user site-packages directory:

    CXC_HOOK_SRC=/path/to/chimerax_clipper_parallel_build.py \
        ChimeraX --nogui --silent --exit --script install_parallel_build_hook.py

Set CXC_HOOK_UNINSTALL=1 to remove the hook instead of installing it. The target
directory is whatever the running ChimeraX reports as its toolshed site dir, so
this always matches the ChimeraX it is run under (per-version user dir).
"""

import os
import shutil

MODULE_NAME = "chimerax_clipper_parallel_build.py"
PTH_NAME = "chimerax_clipper_parallel_build.pth"


def _site_dir():
    # `session` is injected into the script namespace by ChimeraX's runscript.
    return session.toolshed._site_dir  # noqa: F821


def main():
    site_dir = _site_dir()
    module_dst = os.path.join(site_dir, MODULE_NAME)
    pth_dst = os.path.join(site_dir, PTH_NAME)

    if os.environ.get("CXC_HOOK_UNINSTALL"):
        removed = []
        for p in (module_dst, pth_dst):
            if os.path.exists(p):
                os.remove(p)
                removed.append(p)
        if removed:
            print("removed parallel-build hook:")
            for p in removed:
                print("  ", p)
        else:
            print("parallel-build hook not present in", site_dir)
        return

    src = os.environ.get("CXC_HOOK_SRC")
    if not src or not os.path.isfile(src):
        print("install_parallel_build_hook: CXC_HOOK_SRC not set or not a file:",
              src)
        return

    os.makedirs(site_dir, exist_ok=True)
    shutil.copyfile(src, module_dst)
    with open(pth_dst, "w") as f:
        f.write("import chimerax_clipper_parallel_build\n")
    print("installed parallel-build hook:")
    print("   module:", module_dst)
    print("   pth:   ", pth_dst)


main()
