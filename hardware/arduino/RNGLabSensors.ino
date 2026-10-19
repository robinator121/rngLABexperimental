
#include <Wire.h>

// RNGLab Sensor Acquisition v0.1
// Arduino Uno R3
// MPU6050 + BMP280/BME280 + VS1838B

#define IR_PIN 2
#define SAMPLE_INTERVAL_MS 50

uint8_t mpuAddr = 0;
uint8_t baroAddr = 0;
uint8_t baroType = 0;

volatile uint32_t irPulses = 0;
volatile uint32_t lastIR = 0;

uint32_t lastSample = 0;
uint32_t previousPulses = 0;

void onIR() {
  irPulses++;
  lastIR = micros();
}

bool readRegs(uint8_t addr, uint8_t reg,
              uint8_t *buf, uint8_t count) {
  Wire.beginTransmission(addr);
  Wire.write(reg);

  if (Wire.endTransmission(false) != 0)
    return false;

  if (Wire.requestFrom(addr, count) != count)
    return false;

  for (uint8_t i = 0; i < count; i++)
    buf[i] = Wire.read();

  return true;
}

bool readReg(uint8_t addr, uint8_t reg,
             uint8_t &value) {
  return readRegs(addr, reg, &value, 1);
}

bool writeReg(uint8_t addr, uint8_t reg,
              uint8_t value) {
  Wire.beginTransmission(addr);
  Wire.write(reg);
  Wire.write(value);
  return Wire.endTransmission() == 0;
}

void initMPU() {
  uint8_t id;

  for (uint8_t a = 0x68; a <= 0x69; a++) {
    if (readReg(a, 0x75, id) && id == 0x68) {
      mpuAddr = a;
      writeReg(a, 0x6B, 0x00);
      writeReg(a, 0x1B, 0x00);
      writeReg(a, 0x1C, 0x00);
      break;
    }
  }

  Serial.print(F("# MPU6050: "));
  Serial.println(mpuAddr ? F("FOUND") :
                           F("NOT FOUND"));
}

void initBarometer() {
  uint8_t id;

  for (uint8_t a = 0x76; a <= 0x77; a++) {
    if (!readReg(a, 0xD0, id))
      continue;

    if (id == 0x58 || id == 0x60 ||
        id == 0x56 || id == 0x57) {
      baroAddr = a;
      baroType = id;
      break;
    }
  }

  if (!baroAddr) {
    Serial.println(F("# BAROMETER: NOT FOUND"));
    return;
  }

  Serial.print(F("# BAROMETER: "));
  Serial.println(baroType == 0x60 ?
                 F("BME280") : F("BMP280"));

  // Disable internal IIR filtering.
  writeReg(baroAddr, 0xF5, 0x00);

  // Enable humidity sampling on BME280.
  if (baroType == 0x60)
    writeReg(baroAddr, 0xF2, 0x01);

  // Temperature x1, pressure x1, normal mode.
  writeReg(baroAddr, 0xF4, 0x27);
}

void printMPU() {
  uint8_t b[14];

  if (!mpuAddr || !readRegs(mpuAddr, 0x3B, b, 14)) {
    Serial.print(F("NA,NA,NA,NA,NA,NA,"));
    return;
  }

  const uint8_t offsets[] = {0, 2, 4, 8, 10, 12};

  for (uint8_t i = 0; i < 6; i++) {
    uint8_t k = offsets[i];

    int16_t value = (int16_t)(
      ((uint16_t)b[k] << 8) | b[k + 1]
    );

    Serial.print(value);
    Serial.print(',');
  }
}

void printBarometer() {
  uint8_t b[8];
  uint8_t count = baroType == 0x60 ? 8 : 6;

  if (!baroAddr ||
      !readRegs(baroAddr, 0xF7, b, count)) {
    Serial.print(F("NA,NA,NA,"));
    return;
  }

  uint32_t pressure =
    ((uint32_t)b[0] << 12) |
    ((uint32_t)b[1] << 4) |
    (b[2] >> 4);

  uint32_t temperature =
    ((uint32_t)b[3] << 12) |
    ((uint32_t)b[4] << 4) |
    (b[5] >> 4);

  Serial.print(pressure);
  Serial.print(',');
  Serial.print(temperature);
  Serial.print(',');

  if (baroType == 0x60) {
    uint16_t humidity =
      ((uint16_t)b[6] << 8) | b[7];
    Serial.print(humidity);
  } else {
    Serial.print(F("NA"));
  }

  Serial.print(',');
}

void printIR() {
  uint32_t pulses;
  uint32_t last;

  noInterrupts();
  pulses = irPulses;
  last = lastIR;
  interrupts();

  Serial.print(pulses - previousPulses);
  Serial.print(',');

  previousPulses = pulses;

  if (pulses > 0) {
    Serial.println((uint32_t)(micros() - last));
  } else {
    Serial.println(F("NA"));
  }
}

void setup() {
  Serial.begin(115200);
  Wire.begin();
  Wire.setClock(100000);

  delay(500);

  Serial.println(F("# RNGLab v0.1 START"));

  initMPU();
  initBarometer();

  pinMode(IR_PIN, INPUT_PULLUP);
  attachInterrupt(
    digitalPinToInterrupt(IR_PIN),
    onIR,
    FALLING
  );

  Serial.println(
    F("time_us,ax,ay,az,gx,gy,gz,"
      "pressure_raw,temp_raw,humidity_raw,"
      "ir_pulses,ir_age_us")
  );
}

void loop() {
  uint32_t now = millis();

  if ((uint32_t)(now - lastSample) <
      SAMPLE_INTERVAL_MS)
    return;

  lastSample = now;

  Serial.print(micros());
  Serial.print(',');

  printMPU();
  printBarometer();
  printIR();
}
