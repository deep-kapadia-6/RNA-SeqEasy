from app.config import create_default_config, save_config_to_yaml, load_config_from_yaml

# Create default config
config = create_default_config(
    input_file="test_data/example.h5ad",
    output_dir="outputs/test",
    metadata_file="test_data/metadata.csv",
)

# Save to YAML
save_config_to_yaml(config, "configs/test_config.yaml")

# Load from YAML
loaded_config = load_config_from_yaml("configs/test_config.yaml")

print("Config created and loaded successfully!")
print(f"Input file: {loaded_config.input.input_file}")
print(f"QC min genes: {loaded_config.qc.min_genes_per_cell}")
