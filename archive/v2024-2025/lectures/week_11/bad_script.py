def calculate_average(numbers):
    total = sum(numbers)
    count = len(numbers)
    average = total / count
    return average

# Deliberate error: This list has a string that will cause a TypeError
data = [10, 20, 30, "40", 50]
result = calculate_average(data)
print(f"The average is {result}")
