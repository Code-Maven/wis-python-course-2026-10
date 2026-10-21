import csv
from pathlib import Path


def calculate_bmi(height_cm: float, weight_kg: float) -> float:
    """Calculate BMI given height in centimeters and weight in kilograms."""
    height_m = height_cm / 100.0
    return round(weight_kg / (height_m ** 2), 2)


def categorize_bmi(bmi: float) -> str:
    """Categorize a BMI value into standard health categories."""
    if bmi < 18.5:
        return "Underweight"
    elif bmi < 25.0:
        return "Healthy"
    elif bmi < 30.0:
        return "Overweight"
    else:
        return "Obesity"


def process_people_csv(
    input_path: Path = Path("people.csv"),
    output_path: Path = Path("people_with_bmi.csv"),
) -> list[dict[str, str | float]]:
    rows: list[dict[str, str | float]] = []

    with input_path.open(mode="r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if not row.get("ID"):
                continue
            height_cm = float(row["Height"])
            weight_kg = float(row["Weight"])
            bmi = calculate_bmi(height_cm, weight_kg)
            category = categorize_bmi(bmi)

            row["BMI"] = bmi
            row["Category"] = category
            rows.append(row)

    fieldnames = ["ID", "Height", "Weight", "BMI", "Category"]
    with output_path.open(mode="w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main() -> None:
    input_file = Path(__file__).parent / "people.csv"
    output_file = Path(__file__).parent / "people_with_bmi.csv"

    rows = process_people_csv(input_file, output_file)

    print(f"{'ID':<4} {'Height':<8} {'Weight':<8} {'BMI':<8} {'Category'}")
    print("-" * 42)
    for row in rows:
        print(
            f"{row['ID']:<4} {row['Height']:<8} {row['Weight']:<8} "
            f"{row['BMI']:<8.2f} {row['Category']}"
        )
    print(f"\nSaved updated data to {output_file}")


if __name__ == "__main__":
    main()
