#include <ESP8266WiFi.h>
#include <ESP8266HTTPClient.h>
#include "xOD01.h"

const char* ssid = "MNet-guest";               // your network SSID (name)
const char* pass = "LionsInKenya";             // your network password
char servername[]="celestrak.com";           // Celestrak Server
int TLE_select = 6;

xOD01 OD01;

const int DELAY_TIME = 200;

WiFiClient client;

void setup() {
// Starts the I2C communication
#ifdef ESP8266
  Wire.pins(2, 14);
#endif
  Wire.begin();
  
  // Start the OLED Display OD01
  OD01.begin();

  // Small delay
  delay(DELAY_TIME);
  
  Serial.begin(115200);
  Serial.println("Attempting to connect to WiFi");
  OD01.println("Attempting to connect to WiFi");
  WiFi.begin(ssid, pass);
  while ( WiFi.status() != WL_CONNECTED) {
    delay(1000);
    Serial.println("...");
  }

    Serial.println("Connected to wifi");
    OD01.println("Connected to wifi");
    Serial.println("\nStarting connection with server...");
    OD01.println("Starting connection with server...");
  
  makeRequest();
}

void scrollDown(){
  int lineCounter=0;
  char c;
  if (TLE_select==1){
    return;
  }
   while (true) {
   while (!client.available()){
      // while loop runs while waiting for server availability
    }
    c = client.read();
    
    if (c == '\n'){
      lineCounter = lineCounter+1;
    }

    if (lineCounter==3*(TLE_select-1)){
      break;
    }
  }
  return;
}

void makeRequest(){
    // if you get a connection, report back via serial:
    if (client.connect(servername, 80)) {
    Serial.println("connected to server");
    OD01.println("connected to server");
    Serial.println();
    Serial.print("TLE for: ");
    // Make HTTP request:
    client.println("GET /NORAD/elements/stations.txt HTTP/1.0");     // rest of url for your chosen txt file, i.e extension following celestrak.com , Replace everything EXCEPT: GET HTTP/1.0
    client.println();                                                         
    }
    
   // if there are incoming bytes available
   // from the server, read them and print them:
  char c;
  int lineCounter=0; 
 while (!client.available()){
  // while loop runs while waiting for server availability
 }

// Skip HTTP headers
 char endOfHeaders[] = "\r\n\r\n";
  if (!client.find(endOfHeaders))
  {
    Serial.println(F("Invalid response"));
    OD01.println(F("Invalid response"));
    return;
  }
scrollDown();
 lineCounter=0; 
 while (true) {
   while (!client.available()){
   // while loop runs while waiting for server availability
    }
    c = client.read();
    Serial.print(c);
    
    if (c == '\n'){
      lineCounter = lineCounter+1;
    }

    if (lineCounter==3){
      client.stop();
      break;
    }
  }

  // if the server becomes disconnected, stop the client:
  if (!client.connected()) {
    Serial.println();
    Serial.println("disconnecting from server");
    OD01.println("disconnecting from server");
    client.stop();
  }
}

void loop() {
}
