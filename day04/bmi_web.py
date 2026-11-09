from flask import Flask, render_template_string, request

from bmi import calculate_bmi

app = Flask(__name__)

FORM_TEMPLATE = """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>BMI Calculator</title>
</head>
<body>
  <h1>BMI Calculator</h1>
  <form method="post">
    <label for="height_cm">Height (cm):</label>
    <input id="height_cm" name="height_cm" type="number" step="0.01" required>
    <br>
    <label for="weight_kg">Weight (kg):</label>
    <input id="weight_kg" name="weight_kg" type="number" step="0.01" required>
    <br>
    <button type="submit">Calculate BMI</button>
  </form>

  {% if error %}
    <p style="color: red;">{{ error }}</p>
  {% endif %}

  {% if bmi is not none %}
    <p>Your BMI is: <strong>{{ bmi }}</strong></p>
  {% endif %}
</body>
</html>
"""


@app.route("/", methods=["GET", "POST"])
def index() -> str:
    bmi = None
    error = None

    if request.method == "POST":
        try:
            height_cm = float(request.form["height_cm"])
            weight_kg = float(request.form["weight_kg"])
            if height_cm <= 0 or weight_kg <= 0:
                raise ValueError
            bmi = calculate_bmi(height_cm, weight_kg)
        except (KeyError, ValueError):
            error = "Please enter valid positive numbers for height and weight."

    return render_template_string(FORM_TEMPLATE, bmi=bmi, error=error)


if __name__ == "__main__":
    app.run(debug=True)
