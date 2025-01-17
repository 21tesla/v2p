import pandas as pd
import os, sys

expected = pd.read_csv(os.getenv('V2P_DIR') + '/test/expected.csv')
actual = pd.read_csv(sys.argv[1])

assert(expected.equals(actual))
print("Generated predictions matched the expected predictions")
