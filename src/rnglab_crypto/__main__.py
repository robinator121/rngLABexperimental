"""CPU-friendly CLI; no network targets, credentials or real-key inputs."""
import argparse
import json
from pathlib import Path

from rnglab_crypto.experiments.config import CampaignConfig
from rnglab_crypto.experiments.campaign import run_campaign
from rnglab_crypto.reporting import freeze_run, finish_run, verify, write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run")
    run.add_argument("--output", type=Path, required=True)
    choice = run.add_mutually_exclusive_group(required=True)
    choice.add_argument("--config", type=Path)
    choice.add_argument("--smoke", action="store_true")
    run.add_argument("--no-plots", action="store_true")
    check = commands.add_parser("verify")
    check.add_argument("output", type=Path)
    physical = commands.add_parser("physical-audit")
    physical.add_argument("input", type=Path)
    physical.add_argument("--metadata", type=Path, required=True)
    physical.add_argument("--output", type=Path, required=True)
    os_ref = commands.add_parser("os-reference")
    os_ref.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "verify":
        print(json.dumps(verify(args.output), sort_keys=True))
    elif args.command == "physical-audit":
        from rnglab_crypto.sources.physical import audit_waveforms
        if args.output.exists():
            parser.error("output exists; choose a new audit path")
        data = audit_waveforms(args.input, json.loads(args.metadata.read_text()))
        write_json(args.output, data)
    elif args.command == "os-reference":
        from rnglab_crypto.sources import OSReference
        from rnglab_crypto.monitoring import diagnostics
        from rnglab_crypto.reporting import digest, environment
        args.output.mkdir(parents=True, exist_ok=False)
        source = OSReference()
        data = source.read_bytes(4096)
        (args.output/"public-sample.bin").write_bytes(data)
        write_json(args.output/"reference.json", {"provenance": source.manifest.to_dict(),
                   "sha256": digest(args.output/"public-sample.bin"), "diagnostics": diagnostics(data),
                   "environment": environment(), "key_generation": False,
                   "qualification": "one observational OS sample; public after capture; not certification"})
    else:
        config = CampaignConfig.smoke() if args.smoke else CampaignConfig(**json.loads(args.config.read_text()))
        out = freeze_run(args.output, config)
        results, performance = run_campaign(config, lambda msg: print(msg, flush=True))
        finish_run(out, results, performance, make_plots=not args.no_plots)
        print(json.dumps(verify(out), sort_keys=True))


if __name__ == "__main__":
    main()
