import argparse
import json
from pathlib import Path
from rnglab.experiments.runner import run


def main():
    parser = argparse.ArgumentParser(description="Reproducible RNG degradation / Monte Carlo prototype")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    try:
        config = json.loads(args.config.read_text())
        if args.output.exists() and any(args.output.iterdir()):
            parser.error("output directory must be empty (preserves existing results)")
        data = run(config, args.output)
        if not args.no_plots:
            from rnglab.experiments.plotting import plot_all
            plot_all(config, args.output, *data)
    except (ValueError, KeyError) as exc:
        parser.error(str(exc))
    print(f"Results saved to {args.output.resolve()}")


if __name__ == "__main__":
    main()
