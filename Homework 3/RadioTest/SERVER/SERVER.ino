/*************************************************************
  This is an examples for the RL01-02-03 Radio Range

  You can buy one on our store!
  -----> https://xinabox.cc/products/RL01/
  -----> https://xinabox.cc/products/RL02/
  -----> https://xinabox.cc/products/RL03/

  This example requests the Alcohol sensor to measure
  the Breath Alcohol Level

  Currently Supported on the following ☒CHIPs:
  - CW01
  - CR01/02/03

  The sensor communicates over the I2C Bus.

*************************************************************/

#include <xCore.h>
#include <xRL0x.h>

#define RL03_FREQ 915.0

#define CW01_RED 12
#define CW01_GREEN 13
#define CW01_BLUE 5

int cycle=0;
int count=0;
int flag=0;

void setup() {
  // Start the Serial Monitor
  Serial.begin(115200);

  // Set the RGB Pin directions
  pinMode(LED_RED, OUTPUT);
  pinMode(LED_GREEN, OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);

  // Start the I2C Comunication
  Wire.begin();

  if (!RL0X.begin()) { // <-- enter radio name here
    Serial.println("Check the connector to RL03");
    while (1) {
      // Flash RED to indicate failure
      digitalWrite(LED_RED, HIGH);
      delay(100);
      digitalWrite(LED_RED, LOW);
      delay(100);
    }
  } else {
    // RL0X Initialized correctly
    RL0X.setModemConfig(RL0X.Bw31_25Cr48Sf512);
    RL0X.setFrequency(RL03_FREQ);
    RL0X.setTxPower(23, false);
  }
  Serial.println("setup done");
}

void loop() {
  digitalWrite(LED_BUILTIN,HIGH);
  Serial.println("Waiting");
  if (cycle%4==0){
    RL0X.setModemConfig(RL0X.Bw31_25Cr48Sf512);
  }
  if (cycle%4==1){
    RL0X.setModemConfig(RL0X.Bw500Cr45Sf128);
  }
  if (cycle%4==2){
    RL0X.setModemConfig(RL0X.Bw31_25Cr48Sf512);
  }
  if (cycle%4==3){
    RL0X.setModemConfig(RL0X.Bw500Cr45Sf128);
  }
  if (RL0X.waitAvailableTimeout(6000)) {
    uint8_t buf[195];
    uint8_t len = sizeof(buf);
    if (RL0X.recv(buf, &len)) {
      count=count+1;
      flag=1;
      digitalWrite(LED_RED, HIGH);
      Serial.println("got message: ");
      Serial.println((char*)buf);
      Serial.print("RSSI: ");
      Serial.println(RL0X.lastRssi(), DEC);

      // Send a reply
      uint8_t data[] = "And hello back to you";
      delay(100);
      RL0X.send(data, sizeof(data));
      Serial.println("Sent a reply");
      digitalWrite(LED_RED, LOW);
    } else {
      Serial.println("recv failed");
    }
  }
  digitalWrite(LED_BUILTIN,LOW);
  if (count%10==0 && flag==1){
    cycle=cycle+1;
    Serial.println("Switching Modem Config");
    flag=0;
  }
}
