import re

input_file = "proposal_targets.txt"   # the big block you pasted
output_file = "tic_ids_unique.txt"

tic_ids = set()

with open(input_file, "r") as f:
    for line in f:
        if line.strip().startswith("#"):  # skip comments
            continue
        m = re.match(r"^(\d+)", line.strip())
        if m:
            tic_ids.add(int(m.group(1)))

with open(output_file, "w") as f:
    for tid in sorted(tic_ids):
        f.write(f"{tid}\n")

print(f"Found {len(tic_ids)} unique TIC IDs. Saved to {output_file}")