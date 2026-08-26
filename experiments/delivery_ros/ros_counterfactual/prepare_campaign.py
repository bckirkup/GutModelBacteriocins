import json
from pathlib import Path

ROOT = Path("/home/ubuntu/gutibm-campaign/ros-counterfactual")
CEILING_BASE = Path(
    "/home/ubuntu/gutibm-campaign/ceiling-diag/ceiling_res2/input.json"
)
FLORA_BASE = Path(
    "/home/ubuntu/gutibm-campaign/ceiling-diag/flora_res2/input.json"
)


def prepare(base_path, name, grid_dx, vbf_sink, count, ros):
    config = json.loads(base_path.read_text())
    config["total_time"] = 21600
    config["output_interval"] = 300
    config["grid_dx"] = grid_dx
    config["oxygen.vbf_sink"] = vbf_sink
    config["oxygen.k_ROS"] = ros
    config["hdf5"]["schedule"]["summary"] = 1
    config["hdf5"]["schedule"]["agents"] = 5
    config["hdf5"]["schedule"]["grid"] = 10
    config["hdf5_file"] = "output.h5"
    if count == 4:
        config["initial_strains"][0]["count"] = 4
        config["initial_strains"] = [config["initial_strains"][0]]
    elif count == 80:
        config["initial_strains"][0]["count"] = 60
        config["initial_strains"][1]["count"] = 20
    else:
        raise ValueError(f"unsupported founder count: {count}")
    configured_total = sum(
        strain["count"] for strain in config["initial_strains"]
    )
    if configured_total != count:
        raise ValueError(
            f"{name}: configured founders {configured_total} != {count}"
        )
    output_dir = ROOT / name
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "input.json").write_text(
        json.dumps(config, indent=2) + "\n"
    )


def main():
    prepare(CEILING_BASE, "A_ros0_res2", 2e-6, 0.0, 4, 0.0)
    prepare(CEILING_BASE, "A_ctrl_res2", 2e-6, 0.0, 4, 1e2)
    prepare(FLORA_BASE, "B_ros0_res2", 2e-6, 1e-3, 80, 0.0)
    prepare(FLORA_BASE, "B_ctrl_res2", 2e-6, 1e-3, 80, 1e2)
    prepare(FLORA_BASE, "B_ros0_res6", 6e-6, 1e-3, 80, 0.0)
    prepare(FLORA_BASE, "B_ctrl_res6", 6e-6, 1e-3, 80, 1e2)


if __name__ == "__main__":
    main()
