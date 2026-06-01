import rcpgenerator
from pathlib import Path


def main() -> None:
    packing = rcpgenerator.Packing(
        phi=0.11,
        N=4,
        Ndim=2,
        box=[1.0, 1.0],
        walls=[0, 0],
        fix_height=False,
        dist={"type": "mono", "d": 1.0},
        neighbor_max=0,
        seed=123,
    )

    print("needs_initialize before pack:", packing.needs_initialize)
    print("needs_pack before pack:", packing.needs_pack)
    result = packing.pack()
    image_path = packing.savefig(
        path=Path(__file__).resolve().parent / "output" / "starter.png",
        palette_choice=7,
    )
    print("final phi:", result["phi"])
    print("force samples:", len(result["force_history"]))
    print("image:", image_path)


if __name__ == "__main__":
    main()
