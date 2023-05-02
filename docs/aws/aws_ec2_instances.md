---  
hide:  
  - navigation  
---
# AWS Login

### Connect AWS instance on ubuntu/Mac

  1. Download key 
    - [Linux/Mac OS](https://course-cg-5534.s3.amazonaws.com/aws_key/course-setup-student-key.pem){:target="_blank"}
    - [Windows OS](https://course-cg-5534.s3.amazonaws.com/aws_key/course-setup-student-key.ppk){:target="_blank"}

  2. Locate the private key and set the permissions
    ```
        chmod 400 course-setup-student-key.pem
    ```
  3. Connect to your Linux instance using an SSH client
    ```
        ssh -i course-setup-student-key.pem <Instance ID>
    ```
    <!-- **Note**: `<Instance ID>` - check in the below table  -->


**Note**: 

- Fetch `<Instance ID>` from your name in the below table.
- :fontawesome-regular-clock: **AWS Instance** :<span class="add-info"> Monday - Friday ( 13:00 - 22:00 ) </span>

----

## List of user and instance details
  

|Name|Instance ID|Status|
| :--- | :--- | :--- |
|Caroline Hesselager|-|STOPPED|
|Sulaf Abd Own|-|STOPPED|
|Tina Becirovic|-|STOPPED|
|Anna Thulin|-|STOPPED|
|Chao Zheng|-|STOPPED|
|Lina Stepanauskaite|-|STOPPED|
|Yelin Zhao|-|STOPPED|
|Ziyan Ma|-|STOPPED|
|Tove Ekdahl Hjelm|-|STOPPED|
|Petar Mitev|-|STOPPED|
|Linn Hjelmgren|-|STOPPED|
|Andreas Ekholm|-|STOPPED|
|Donal Barrett|-|STOPPED|
|Blaz Oder|-|STOPPED|
|Ye Yuan|-|STOPPED|
|Elin Barnekow|-|STOPPED|
|Evangelos Doukoumopoulos|-|STOPPED|
|Alen Lovric|-|STOPPED|
|Anna Plym|-|STOPPED|
|Ioannis Zerdes|-|STOPPED|