#include <xCore.h>
#include <xSW01.h>
#include <xSI01.h>
#include <xSL01.h>
#include <xSN01.h>

xSW01 SW01;
xSI01 SI01;
xSL01 SL01;
xSN01 SN01;

int count=0;

const int DELAY_TIME = 1000;

void setup() {
  // Start the Serial Monitor
  Serial.begin(115200);

  // Set the I2C Pins for CW01
  #ifdef ESP8266
    Wire.pins(2, 14);
    Wire.setClockStretchLimit(15000);
  #endif

  // Start the I2C Comunication
  Wire.begin();
  
  // Start the sensors
  SW01.begin();
  SN01.begin();
  SL01.begin();
  SI01.begin();
  
  //Delay for sensor to normalise
  delay(15000);

}

void loop() {

  if (count<12) {

  // Create variables to store the data read from SW01
  float alt;
  alt = 0;
  float humidity;
  humidity = 0;
  float pressure;
  pressure = 0;
  float tempC;
  float tempF;
  tempC = tempF = 0;

  // Read and calculate data from SW01 sensor
  SW01.poll();
  
  // Request SW01 to get the altitude measurement and store in
  // the variables
  alt = SW01.getAltitude(101325);
  humidity = SW01.getHumidity();
  pressure = SW01.getPressure();
  tempC = SW01.getTempC(); // Temperature in Celcuis
  tempF = SW01.getTempF(); // Temperature in Farenheit

  if (count>10) {
  Serial.print("SW01: Altitude: ");
  Serial.print(alt);
  Serial.print(" m    ");
  Serial.print("Humidity: ");
  Serial.print(humidity);
  Serial.print("%    ");
  Serial.print("Pressure: ");
  Serial.print(pressure);
  Serial.print(" Pa    ");
  Serial.print("Temperature: ");
  Serial.print(tempC);
  Serial.print(" C    ");
  Serial.print("Temperature: ");
  Serial.print(tempF);
  Serial.println(" F    ");
  }

  // Read and calculate data from SL01 sensor
  SI01.poll();

  if (count>10) {
  Serial.print("SI01: ");
  printGyro();  // Print "G: gx, gy, gz"
  printAccel(); // Print "A: ax, ay, az"
  printMag();   // Print "M: mx, my, mz"
  printAttitude(); // Print Roll, Pitch and G-Force
  }

  // Create a variables to store the incoming measurements
  // from SL01 sensor
  float lux;
  lux = 0;
  float uv;
  uv = 0;
  
  // Poll Sensor for collect data
  SL01.poll();

  // Request SL01 to return calculated LUX intensity
  lux = SL01.getLUX();

  if (count>10) {
  // Display Data on the Serial monitor
  Serial.print("SL01: Ambient Light Level: ");
  Serial.print(lux);
  Serial.print(" LUX    ");

  uv = SL01.getUVA();
  Serial.print("UVA Intersity: ");
  Serial.print(uv);
  Serial.print(" uW/m^2    ");
  uv = SL01.getUVB();
  Serial.print("UVB Intensity: ");
  Serial.print(uv);
  Serial.print(" uW/m^2    ");
  uv = SL01.getUVIndex();
  Serial.print("UVB Index: ");
  Serial.println(uv);
  }

  // Create a variable to store the data read from SN01
  String time;
  long latitude = 0;
  long longitude = 0;
  String date;

  // Poll the sensor to read all available data
  SN01.poll();


  // Get the date from GPS
  date = SN01.getDate();
  // Get the time from the GPS 
  time = SN01.getTime();
  // Get the latitude from GPS
  latitude = SN01.getLatitude();
  // Get the longitude from GPS
  longitude = SN01.getLongitude();

  if (count>10) {
  // Display the recorded data over the serial monitor
  Serial.print("SN01: GPS Time: ");
  Serial.print(time);
  Serial.print("    ");
  Serial.print("GPS Date: ");
  Serial.print(date);
  Serial.print("    ");
  Serial.print("GPS Latitude: ");
  Serial.print(latitude);
  Serial.print("    ");
  Serial.print("GPS longitude: ");
  Serial.println(longitude);

  Serial.println();
  }
  }

  // Small delay between sensor reads
  delay(DELAY_TIME);
  count=count+1;
  
}

void printGyro(void) {
  Serial.print("G: [");
  Serial.print(SI01.getGX(), 2);
  Serial.print(", ");
  Serial.print(SI01.getGY(), 2);
  Serial.print(", ");
  Serial.print(SI01.getGZ(), 2);
  Serial.print("]    ");

}

void printAccel(void) {
  Serial.print("A: [");
  Serial.print(SI01.getAX(), 2);
  Serial.print(", ");
  Serial.print(SI01.getAY(), 2);
  Serial.print(", ");
  Serial.print(SI01.getAZ(), 2);
  Serial.print("]    ");
}

void printMag(void) {
  Serial.print("M: [");
  Serial.print(SI01.getMX(), 2);
  Serial.print(", ");
  Serial.print(SI01.getMY(), 2);
  Serial.print(", ");
  Serial.print(SI01.getMZ(), 2);
  Serial.print("]    ");

}

void printAttitude(void) {
  Serial.print("Roll: ");
  Serial.print(SI01.getRoll(), 2);
  Serial.print("    ");
  Serial.print("Pitch :");
  Serial.print(SI01.getPitch(), 2);
  Serial.print("    ");
  Serial.print("GForce :");
  Serial.print(SI01.getGForce(), 2);
  Serial.println("    ");
}
