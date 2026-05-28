# ================================================================================

# Google Maps API routing to capture travel distance and time for the shortest path
# Input data:  id	respondent	latO	lgnO	latD	lgnD	category
# Output data: respondent	rid	o_lat	o_lon	category	d_lat	d_lon	distance_m	travel_time_min

# ================================================================================

import os
import time
import requests
import pandas as pd
from pathlib import Path

# ==============================
# CONFIGURATION
# ==============================

# Set your Google Maps API key here

API_KEY = 'xxx'     # Provide your API key

BASE_DIR = Path(r"E:\xxx")      # Add path to data

ORIGIN_FILE = BASE_DIR / "origin_destination.csv"     # Update database names
OUTPUT_FILE = BASE_DIR / "impedance_mode_D_T.csv"     # Update name of output database     

# Pause between API calls to avoid rate limits
SLEEP_BETWEEN_CALLS = 0.05


# ==============================
# GOOGLE MAPS - ROUTING 
# ==============================

def get_distance_and_time(origin_lat, origin_lng, dest_lat, dest_lng):
    url = "https://maps.googleapis.com/maps/api/distancematrix/json"
    params = {
        "origins": f"{origin_lat},{origin_lng}",
        "destinations": f"{dest_lat},{dest_lng}",
        "mode": "driving",
        "units": "metric",
        "key": API_KEY,
    }

    response = requests.get(url, params=params)
    data = response.json()

    if data.get("status") != "OK":
        raise RuntimeError(f"Top-level API error: {data.get('status')}")

    element = data["rows"][0]["elements"][0]

    # IMPORTANT: return if OK
    if element.get("status") != "OK":
        return None, None

    distance_m = element["distance"]["value"]
    duration_s = element["duration"]["value"]
    duration_min = round(duration_s / 60, 2)

    return distance_m, duration_min


# ==============================
# LOOP FOR ROWS
# ==============================

def main():
    data = pd.read_csv(ORIGIN_FILE, sep=';')
    results = []

    for _, data_row in data.iterrows():
        rid = data_row["id"]
        category = data_row["category"]
        respondent = data_row["respondent"]
        o_lat = data_row["latO"]
        o_lon = data_row["lgnO"]
        d_lat = data_row["latD"]
        d_lon = data_row["lgnD"]

        distance_m = None
        duration_min = None

        try:
            distance_m, duration_min = get_distance_and_time(o_lat, o_lon, d_lat, d_lon)
        except Exception as e:
            print(f"Error for respondent {respondent}, rid {rid}: {e}")

        results.append({
            "respondent": respondent,
            "rid": rid,
            "o_lat": o_lat,
            "o_lon": o_lon,
            "category": category,
            "d_lat": d_lat,
            "d_lon": d_lon,
            "distance_m": distance_m,
            "travel_time_min": duration_min
        })

        time.sleep(SLEEP_BETWEEN_CALLS)

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Saved {len(results_df)} rows to {OUTPUT_FILE}")

    return results_df  # optional


if __name__ == "__main__":
    main()

