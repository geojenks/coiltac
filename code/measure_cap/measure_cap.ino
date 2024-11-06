#include <Wire.h>

const int OUT_PIN = A2;
const int IN_PIN = A0;
const float IN_CAP_TO_GND = 15.0426; // to be experimentally validated
const float R_PULLUP = 34.8;
const int MAX_ADC_VALUE = 4095;
float last_cap = 0;
//const float offset = 4 / 1.5; // if needed
// May add more channels like so:
//const int OUT_PIN2 = A1;
//const int IN_PIN2 = A3;
//const float IN_CAP_TO_GND2 = 15.0426; // to be experimentally validated
//const float offset2 = 4 / 1.5; // if needed
//float last_cap2 = 0;
long int t0;

void setup()
{
  analogReadResolution(12); // Set the analog resolution to 12 bits
  pinMode(OUT_PIN, OUTPUT);
  pinMode(IN_PIN, OUTPUT);
  Serial.begin(115200);
  //Serial.println("***RESET***");
  delay(2000);
  t0 = millis();
}

void loop()
{
  pinMode(IN_PIN, INPUT);
  digitalWrite(OUT_PIN, HIGH);
  int val = analogRead(IN_PIN);
  digitalWrite(OUT_PIN, LOW);
  
  // for more channels...
  //pinMode(IN_PIN2, INPUT);
  //digitalWrite(OUT_PIN2, HIGH);
  //int val2 = analogRead(IN_PIN2);
  //Serial.println(val);
  //digitalWrite(OUT_PIN2, LOW);
  
// Threshold based on RC time constant and measurement rate - generally not needed
//  if (val < 4000)
//  {
    pinMode(IN_PIN, OUTPUT);
    float capacitance = (float)val * IN_CAP_TO_GND / (float)(MAX_ADC_VALUE - val);// - offset;
    //Serial.print("Cap = ");
    //Serial.println(val);
    //Serial.print(", ");
    //Serial.print(millis()-t0);
    //Serial.print(", ");
    //exponential smoothing if wanted:
    //if (last_cap != 0){capacitance = 0.99*last_cap + 0.01* capacitance;}
    Serial.println(capacitance, 5);
    last_cap = capacitance;
    // for more channels...
    //pinMode(IN_PIN2, OUTPUT);
    //float capacitance2 = (float)val2 * IN_CAP_TO_GND / (float)(MAX_ADC_VALUE - val);// - offset;
    //Serial.println(capacitance2, 5);
    //Serial.println("pF");
    //last_cap2 = capacitance2;
    ;
//  }
/*
 //  Section only needed for significantly larger capacitance, where caps do not dishcarge quickly enough and RC must 
 // be calculated instead of Cap directly
  else
  {
    pinMode(IN_PIN, OUTPUT);
    delay(1);
    pinMode(OUT_PIN, INPUT_PULLUP);
    unsigned long u1 = micros();
    unsigned long t;
    int digVal;

    Serial.println("pre");
    do
    {
      digVal = digitalRead(OUT_PIN);
      unsigned long u2 = micros();
      t = u2 > u1 ? u2 - u1 : u1 - u2;
    }
    while ((digVal < 1) && (t < 400000L));
    Serial.println("post");

    pinMode(OUT_PIN, INPUT);
    val = analogRead(OUT_PIN);
    digitalWrite(IN_PIN, HIGH);
    int dischargeTime = (int)(t / 1000L) * 5;
    delay(dischargeTime);
    pinMode(OUT_PIN, OUTPUT);
    digitalWrite(OUT_PIN, LOW);
    digitalWrite(IN_PIN, LOW);

    float capacitance = -(float)t / R_PULLUP / log(1.0 - (float)val / (float)MAX_ADC_VALUE);

    if (capacitance > 1000.0)
    {
      Serial.print(capacitance / 1000.0, 2);
      Serial.println("uF");
    }
    else
    {
      Serial.print(capacitance);
      Serial.println("nF");

    }
  last_cap = capacitance;
  }
*/
  //while (micros() % 1000 != 0);
}
